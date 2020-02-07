# Create an executable file: pyinstaller --onefile pathologist.spec
import numpy as np
import pandas as pd
import time
import multiprocessing as mp
import sys
import gc
try:
    from udp import calc_udp_gmm, calc_udp_nbm
except ImportError:
    from .udp import calc_udp_gmm, calc_udp_nbm
from scipy.stats import mannwhitneyu as mann


class path_activity:
    def __init__(self, udp, is_rnaseq=False):
        print(time.ctime(), 'Init activity object')
        self.orig_paths = pd.read_csv('pathologist.db.txt', delimiter='\t', header=None)
        self.orig_paths.columns = ['path_name', 'path_id', 'molType', 'molName', 'molNum', 'molLink', 'c7', 'c8', 'c9',
                                   'molRole', 'intID', 'intType']
        self.orig_paths.drop(['c8', 'c9'], inplace=True, axis=1)
        self.orig_paths = self.orig_paths.sort_values(['path_id', 'intID'])
        up2ll = pd.read_csv('UP2LL.txt', delimiter='\t', header=None, index_col=0)
        self.up2ll_dict = up2ll[1].to_dict()
        self.orig_cmplx = pd.read_csv('pathologist.complexes.txt', delimiter='\t', header=None)
        # Columns: (0==complex ID, 1==molecule type, 3==molecule ID, 4==Links).
        self.orig_cmplx.drop([2, 5, 6, 7], inplace=True, axis=1)
        self.orig_cmplx.columns = ['complex_ID', 'molType', 'molID', 'molLink']
        self.is_rnaseq = is_rnaseq

        if (is_rnaseq):
            # For RNASeq data we need to create a probe to gene dictionary using Affymetrix table.
            self.probe_to_gene_df = pd.read_csv('HG-U133_Plus_2.na36.annot.csv', index_col=0)
            # Split 'Gene Symbol' in case of multiple genes and stack them in a series.
            new_col = self.probe_to_gene_df['GeneSymbol'].str.split('|', expand=True).stack()
            # Join the new series (after setting a name and index) back to the dataframe.
            self.probe_to_gene_df = self.probe_to_gene_df.join(
                pd.Series(index=new_col.index.droplevel(1), data=new_col.values, name='gene'))
            self.probe_to_gene_df.drop_duplicates(inplace=True)
            del self.probe_to_gene_df['GeneSymbol']
            self.probe_to_gene_df.reset_index(inplace=True)
            self.probe_to_gene_df.columns = ['probe', 'gene']

        self.udp = udp
        self.probelinks = pd.read_csv('probelinksv2.txt', index_col=0) #, header=None
        #self.probelinks.drop(2, inplace=True, axis=1)
        #self.probelinks.columns = ['probe', 'link']
        self.probelinks['link'] = self.probelinks.link.replace('---', 0)
        self.probelinks['link'] = pd.to_numeric(self.probelinks.link)
        # Use only probes or genes that appear in paths.
        #self.filter_probes()

    # Create dictionaries mapping links to the relevant probabilities (UDP), for specific sample_num.
    def calc_link_to_pr(self, sample):
        udp_sample = self.udp.loc[:, sample].copy().reset_index()
        udp_sample.columns = ['probe', 'pr']
        probe_to_pr = pd.merge(self.probelinks, udp_sample, how='left')
        probe_to_pr = probe_to_pr.groupby(['link'], sort=False)['pr'].max()
        self.link_to_pr_dict = probe_to_pr.to_dict()

    def calc_link_to_gene_to_pr(self, sample):
        udp_sample = self.udp.loc[:, sample].copy().reset_index()
        udp_sample.columns = ['gene', 'pr']
        # udp_sample['gene'] = udp_sample.gene.apply(lambda x: x.split('|')[0]) #Remove Entre gene IDs.
        #probelinks['gene'] = probelinks.probe.apply(lambda x: probe_to_gene_dict.get(x, None))
        probe_to_pr = pd.merge(
            self.probelinks, self.probe_to_gene_df, on='probe', how='left')
        probe_to_pr = pd.merge(probe_to_pr, udp_sample, how='left')
        probe_to_pr = probe_to_pr.groupby(['link'], sort=False)['pr'].max()
        self.link_to_pr_dict = probe_to_pr.to_dict()

    # Replace link with UDP.
    def link_to_udp(self, link):
        if ('LL' in link):
            return (self.link_to_pr_dict.get(int(link.replace('LL:', '')), np.NaN))
        if ('UP' in link):
            return (self.link_to_pr_dict.get(int(self.up2ll_dict.get(link.replace('UP:', ''), 0)), np.NaN))
        return np.NaN

    # Remove link prefix.
    def rem_link_prefix(self, link):
        if ('LL' in link):
            return (int(link.replace('LL:', '')))
        if ('UP' in link):
            return (int(self.up2ll_dict.get(link.replace('UP:', ''), 0)))
        return np.NaN

    # Filter only probes that appear in paths. We collect all links that are mentioned in either the paths or complexes database.
    def filter_probes(self):
        path_links = set()

        # Find all links that are used in paths.
        proteins = self.orig_paths.loc[self.orig_paths.molType == 'protein'].copy()
        proteins.molLink.apply(lambda x: [path_links.add(self.rem_link_prefix(i)) for i in str(x).split(',')])
        # Find all links that are used in complexes.
        proteins = self.orig_cmplx.loc[self.orig_cmplx.molType == 'protein'].copy()
        proteins.molLink.apply(lambda x: [path_links.add(self.rem_link_prefix(i)) for i in str(x).split(',')])

        path_links.remove(np.nan)
        #new_index = pd.Series(list(path_links))
        probes = self.probelinks.loc[self.probelinks['link'].isin(path_links)].probe.reset_index(drop=True)
        if (not self.is_rnaseq):
            self.udp = self.udp.loc[list(set(self.udp.index) & set(probes.values))]
        else:
            genes = self.probe_to_gene_df.loc[self.probe_to_gene_df.probe.isin(probes)].gene.copy()
            self.udp = self.udp.loc[self.udp.index.isin(genes)]

    # Complexes are proteins == basic complexes (built from links to probes) or group of proteins.
    # Parse complexes in order to build a dictionary for complex->UDP (similar to the links to UDP dictionary).
    # Problem is that complexes might be composed of other complexes. Hence we need another join operation.
    # We call complexes that have only proteins (i.e. no other complexes) 'basics'. They are part of other (non-basic) complexes.
    # First we calculate UDP for basics, by mutiplyi UDPs of corresponding proteins. Then we can calculate
    # the UDP of all complexes by, again, multiplying all UDPs of mulecules (proteins and basics) in a path.
    # Molecules with no links (i.e. basic complexes) will get a default link. In that case they will get the default UDP value of 1.
    def calc_cmplx_to_pr(self):
        cmplx = self.orig_cmplx.copy()
        # Handle molecules (proteins, rnas and compouns). Parse a link string in the links column in the pathways file or the complexes file.
        cmplx['pr'] = cmplx['molLink'].apply(lambda x: max([self.link_to_udp(i) for i in str(x).split(',')]))
        # Some rows have 'rna' type with no links. They are not really used in complexes.
        cmplx.loc[cmplx['molType'] == 'rna', 'pr'] = 1
        # Compounds are assumed to always be present (UDP == 1).cmplx.loc[cmplx['molType'] == 'compound', 'pr'] = 1

        # Handle basic complexes. The product of molecules on a path are now the accurate basic complex UDP.
        # non-basic complexes don't have links.
        cmplx_basic = cmplx.copy().dropna(subset=['molLink'])
        cmplx_basic = cmplx_basic[['complex_ID', 'pr']].groupby('complex_ID', as_index=False).prod()
        cmplx_basic = cmplx_basic.set_index('complex_ID')
        basic_to_dict = cmplx_basic['pr'].to_dict()

        # Handle non-basic complexes (they can be found by searching a molecule of type 'complex').
        # To handle non-basic complex we just need to mutiply the pr for their group of molecules.
        # Update complx table with the basic complexes values, only for molecules of type 'complex' (this will ignore basic-complex since their column 3 has proteins).
        cmplx.loc[cmplx['molType'] == 'complex', 'pr'] = cmplx['molID'].apply(lambda x: basic_to_dict.get(x, np.NaN))
        # Remove molecules that we couldn't calculate pr for (missing links, unknown basic complexes etc.).
        cmplx = cmplx.dropna()
        cmplx = cmplx[['complex_ID', 'pr']].groupby('complex_ID', as_index=False).prod()
        cmplx = cmplx.set_index('complex_ID')
        return cmplx['pr'].to_dict()

    # Calculate activity and consistency of paths.
    def process_samples(self, udp_chunk):
        # Table columns are (sampleID, path_id, activity, consistency)
        sample_results = np.empty((0, 6))
        for sample in udp_chunk:
            # Calculate UDP of the molecules in all the paths.
            if(not self.is_rnaseq):
                self.calc_link_to_pr(sample)
            else:
                self.calc_link_to_gene_to_pr(sample)
            cmplx_to_pr_dict = self.calc_cmplx_to_pr()
            paths = self.orig_paths.copy()
            # gc.collect()
            # If a molecule does not have a probability we need to remove the whole interaction from activity and consistency calculations.
            paths.loc[paths.molType == 'protein', 'pr'] = paths.molLink.apply(lambda x: max([self.link_to_udp(
                i) for i in str(x).split(',')]))  # Max returns NaN if there is at least one NaN in the list.
            # Compounds are assumed to always be present
            paths.loc[paths.molType == 'compound', 'pr'] = 1
            paths.loc[paths.molType == 'rna', 'pr'] = 1
            paths.loc[paths.molType == 'complex', 'pr'] = paths.molNum.apply(lambda x: cmplx_to_pr_dict.get(x, np.NaN))

            # Calculate activity and consistency of each interaction.
            # Activity of interactions is a multiplication of the inputs.
            paths.loc[paths.molRole == 'inhibitor', 'pr'] = 1 - paths.pr
            # max(paths['path_id'].groupby(paths['intID']).nunique()) #This shows that some interactions participate in few paths.
            # Save the pr for output molecules, since it also receieved a UDP.
            paths['pr_output'] = 1
            paths.loc[paths.molRole == 'output', 'pr_output'] = paths.pr
            paths.loc[paths.molRole == 'output', 'pr'] = 1
            # Handle molecules that we could not find UDP for, we need to remove the whole interaction.
            # Pandas like R is ignoring NA values in gropuby.transform, so instead we zero it.
            #paths = paths.dropna(subset=['pr'])
            paths = paths.fillna({'pr': 0})
            paths['interaction_activity'] = paths.groupby(['path_id', 'intID'])[
                'pr'].transform('prod')
            paths['interaction_consistency'] = paths.interaction_activity * \
                paths.pr_output + \
                (1 - paths.interaction_activity) * (1 - paths.pr_output)

            # Calculate activity and consistency of each path.
            # Both activity and consistency are averages of all interactions activities and consistencies.
            paths['activity'] = paths.groupby(
                'path_id')['interaction_activity'].transform('mean')
            paths['consistency'] = paths.groupby(
                'path_id')['interaction_consistency'].transform('mean')
            paths['sampleID'] = sample
            paths_result = paths[['sampleID', 'path_name', 'path_id',
                                  'activity', 'consistency', 'molRole']].drop_duplicates()
            sample_results = np.concatenate(
                (sample_results, paths_result.values))  # axis=0
        results_df = pd.DataFrame(data=sample_results, columns=[
                                  'sampleID', 'path_name', 'pathID', 'Activity', 'Consistency', 'molRole'])
        return results_df, paths

    # Read chunks of dataframe columns.
    def chunker_columns(self, size):
        len_df = len(self.udp.columns)
        cols = self.udp.columns
        return [cols[pos:min(pos + size, len_df)] for pos in range(0, len_df, size)]

    # Run calc_udp on parallel.
    def calc_activity_consistency_multi_process(self):
        gc.collect()
        print(time.ctime(), 'Calculate activity and consistency...')
        df = pd.DataFrame()
        pool = mp.Pool()  # Use number of CPUs processes.
        results = [pool.apply_async(self.process_samples, args=(x,))
                   for x in self.chunker_columns(20)]
        for p in results:
            df = df.append(p.get()[0])  # f.get(timeout=100)
            print('.', end="")
            sys.stdout.flush()
        df.drop(['molRole'], inplace=True, axis=1)
        df.drop_duplicates(inplace=True)
        df.to_csv('output_activity.csv', index=False)
        self.activity = df
        print(time.ctime(), "Done.")

    def xmlparser(self, path_id, sample_num):
        nodes = set()
        intxns = set()
        colors_green_to_red = ['00FF00', '11FF00', '22FF00', '33FF00', '44FF00', '55FF00', '66FF00', '77FF00', '88FF00',
                               '99FF00', 'AAFF00', 'BBFF00', 'CCFF00', 'DDFF00', 'EEFF00', 'FFFF00', 'FFEE00', 'FFDD00',
                               'FFCC00', 'FFBB00', 'FFAA00', 'FF9900', 'FF8800', 'FF7700', 'FF6600', 'FF5500', 'FF4400',
                               'FF3300', 'FF2200', 'FF1100', 'FF0000']

        print(f'Create Kegg XML for path: {path_id}, sample: {sample_num}')
        if (isinstance(path_id, int)):
            paths = self.orig_paths.loc[self.orig_paths.path_id == path_id].copy()
        else:
            # If instead of path_id the input has the path name.
            paths = self.orig_paths.loc[self.orig_paths.path_name == path_id].copy()
            path_id = paths.iloc[0]['path_id']
        # Add '+' in case of an active molecule.
        paths.loc[paths.c7 == 'active', 'molName'] = paths.molName + '+'

        # We need to get the extra information of interactions activity (that is not saved in the standard calc_activity_consistency run).
        sampleID = self.udp.columns[int(sample_num)]
        _, prbs = self.process_samples([sampleID])
        prbs = prbs.loc[prbs['path_id'] == path_id]
        prbs = prbs[['molNum', 'molRole', 'pr', 'pr_output', 'intID', 'interaction_activity', 'interaction_consistency']]
        paths = pd.merge(paths, prbs, on=['molNum', 'intID', 'molRole'], how='left')

        xml = '<?xml version="1.0"?>' + \
              '<!DOCTYPE pathway SYSTEM "https://www.kegg.jp/kegg/xml/KGML_v0.7.2_.dtd">' + \
              '<pathway name="path:hsa{0}" org="hsa" number="111" title="{1}">'.format(paths.iloc[0, 1], paths.iloc[0, 0])

        for index, row in paths.iterrows():
            # Handle cases of empty molecule name.
            if (isinstance(row['molName'], float)):
                row['molName'] = str(row['molNum'])
            if (row['molName'] not in nodes):
                nodes.add(row['molName'])

                bg_color = "42ECEF"  # Light turquoise
                if (row['molRole'] == 'output'):
                    fg_color = int(np.nan_to_num(row["pr_output"]) * 30)
                    pr = str(np.nan_to_num(row["pr_output"]))
                else:
                    fg_color = int(np.nan_to_num(row["pr"]) * 30)
                    pr = str(np.nan_to_num(row["pr"]))

                if (row['molType'] == 'complex'):
                    # Kegg uses 'group' for complexes.
                    row['molType'] = 'ortholog'
                    bg_color = "C3B6B6"  # Grey
                if (row['molType'] == 'compound'):
                    row['molType'] = 'protein'

                xml += f'<entry id="{row["molNum"]}" name="hsa:4967 hsa:{row["molName"]}" type="{row["molType"]}" ' + \
                       f'link="https://www.kegg.jp/dbget-bin/www_bget?hsa:4967+hsa:{row["molNum"]}"> ' + \
                       f'<graphics name="{row["molName"]}" fgcolor="#{colors_green_to_red[fg_color]}" bgcolor="#{bg_color}" x="{pr}" ' + \
                       f'type="rectangle" width="{min(60, 8 * len(row["molName"]))}" height="17"/> ' + \
                       '</entry>'
            # Add the interaction
            if (row['intID'] not in intxns):
                intxns.add(row['intID'])
                bg_color = "FFF8F8"  # White
                fg_color = int(row["interaction_activity"] * 30)
                pr_act = str(row["interaction_activity"])
                pr_con = str(row["interaction_consistency"])
                # Add the junction representing the interaction.
                xml += f'<entry id="{row["intID"]}" name="{row["intID"]}" type="reaction" ' + \
                       'link="https://www.kegg.jp/dbget-bin/www_bget?C00158"> ' + \
                       f'<graphics name="{row["intType"]}" fgcolor="#{colors_green_to_red[fg_color]}" bgcolor="#{bg_color}" x="{pr_act}" y="{pr_con}" ' + \
                       f'type="circle" width="{min(60, 8 * len(row["intType"]))}" height="15"/> ' + \
                       '</entry>'

            # Add the graph edges.
            if (row['molRole'] == 'input' or row['molRole'] == 'agent'):
                xml += f'<relation entry1="{row["molName"]}" entry2="{row["intID"]}" type="PPrel"> ' \
                    '<subtype name="activation" value="-->"/> ' \
                    '</relation>'
            elif (row['molRole'] == 'output'):
                xml += f'<relation entry1="{row["intID"]}" entry2="{row["molName"]}" type="PPrel"> ' \
                    '<subtype name="activation" value="-->"/> ' \
                    '</relation>'
            elif (row['molRole'] == 'inhibitor'):
                xml += f'<relation entry1="{row["molName"]}" entry2="{row["intID"]}" type="PPrel"> ' \
                    '<subtype name="activation" value="--|"/> ' \
                    '</relation>'
        xml += '</pathway>'
        # _{path_id}_{sample_num}
        text_file = open("output_path.xml", "w")
        text_file.write(xml)
        text_file.close()
        return xml

    def calc_mann_whitney(self, pathname, file1, file2):
        data1 = pd.read_csv(file1)
        data1 = data1.loc[data1.path_name == pathname]
        data2 = pd.read_csv(file2)
        data2 = data2.loc[data2.path_name == pathname]
        n1 = len(data1)
        n2 = len(data2)
        m_u = n1 * n2 * 0.5
        stat, p = mann(data1.Activity, data2.Activity)
        #sigma_u = np.sqrt(n1*n2*(n1+n2+1)*(10/12))
        print(
            f'Stat: {stat}, P-value: {p}, n1: {n1}, n2: {n2}, m_u(expected under H0): {m_u}')


if __name__ == '__main__':
    # activity_obj = path_activity(rma_filename='c:/Work/src/Bayes/Pathweigh/other/firebrowse/ESCA.csv', rnaseq_file=1)  # /data/GSE29013_RMA.csv
    # activity_obj.calc_df_multi_process()
    # activity_obj.calc_activity_consistency_multi_process()
    # activity_obj = path_activity(udp_filename='/data/output_udp.csv', rnaseq_file=0)  # Non RNAseq object
    # activity_obj.xmlparser(1,1)

    udp = pd.read_csv('./data/output_udp.csv', index_col=0)
    activity_obj = path_activity(udp, False)
    # activity_obj.calc_activity_consistency_multi_process()
    activity_obj.xmlparser(100001, 1)

    # Bladder
    #activity_obj.calc_mann_whitney('Bladder cancer(Kegg)', '/data/output_activity_schizophrenia_gse17612.csv', '/data/output_activity_muscle_gse28422.csv')
    # Stat: 2284.0, P-value: 0.029288669514444116, n1: 51, n2: 110, m_u(expected under H0): 2805.0

    #activity_obj.calc_mann_whitney('Bladder cancer(Kegg)', '/data/output_activity_bladder_gse31684.csv', '/data/output_activity_muscle_gse28422.csv')
    # Stat: 2334.0, P-value: 1.3013141984774421e-11, n1: 93, n2: 110, m_u(expected under H0): 5115.0

    #activity_obj.calc_mann_whitney('Bladder cancer(Kegg)', '/data/output_activity_brca_gse50948.csv', '/data/output_activity_muscle_gse28422.csv')
    # Stat: 4124.0, P-value: 2.7848199372065743e-13, n1: 156, n2: 110, m_u(expected under H0): 8580.0

    # BRCA
    #activity_obj.calc_mann_whitney('BRCA pathway reduction(Rotem)', '/data/output_activity_brca_gse50948.csv', '/data/output_activity_muscle_gse28422.csv')
    # Stat: 7658.0, P-value: 0.0679383681844674, n1: 156, n2: 110, m_u(expected under H0): 8580.0
    #activity_obj.calc_mann_whitney('BRCA pathway reduction(Rotem)', '/data/output_activity_brca_gse50948.csv', '/data/output_activity_schizophrenia_gse17612.csv')
    # Stat: 2814.0, P-value: 0.0008647738664741184, n1: 156, n2: 51, m_u(expected under H0): 3978.0

    #activity_obj.calc_mann_whitney('brca1 dependent ub ligase activity(BioCarta)', '/data/output_activity_brca_gse50948.csv', '/data/output_activity_muscle_gse28422.csv')
    # Stat: 8481.0, P-value: 0.4366733613601155, n1: 156, n2: 110, m_u(expected under H0): 8580.0
    # activity_obj.calc_mann_whitney('brca1 dependent ub ligase activity(BioCarta)', '/data/output_activity_brca_gse50948.csv', '/data/output_activity_schizophrenia_gse17612.csv')
    # Stat: 2977.0, P-value: 0.0035279674829234334, n1: 156, n2: 51, m_u(expected under H0): 3978.0

    #activity_obj.calc_mann_whitney('role of brca1 brca2 and atr in cancer susceptibility(BioCarta)', '/data/output_activity_brca_gse50948.csv', '/data/output_activity_muscle_gse28422.csv')
    # Stat: 7729.0, P-value: 0.08434596977542624, n1: 156, n2: 110, m_u(expected under H0): 8580.0
    #activity_obj.calc_mann_whitney('role of brca1 brca2 and atr in cancer susceptibility(BioCarta)', '/data/output_activity_brca_gse50948.csv', '/data/output_activity_schizophrenia_gse17612.csv')
    # Stat: 2109.0, P-value: 2.4322278154806653e-07, n1: 156, n2: 51, m_u(expected under H0): 3978.0
