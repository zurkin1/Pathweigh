#Create an executable file: pyinstaller --onefile pathologist.spec
import numpy as np
import pandas as pd
import time
import sys

probe_to_gene_df = pd.DataFrame()

# Create dictionaries mapping links to the relevant probabilities (UDP), for specific sample_num.
def calc_link_to_pr(udp, sample_num):
    probelinks = pd.read_csv('data/probelinks.txt', delimiter='\t', header=None)
    probelinks.drop(2, inplace=True, axis=1)
    probelinks.columns = ['probe', 'link']
    probelinks['link'] = probelinks.link.replace('---', 0)
    probelinks['link'] = pd.to_numeric(probelinks.link)
    udp_sample = udp.iloc[:, sample_num].copy().reset_index()
    udp_sample.columns = ['probe', 'pr']
    probe_to_pr = pd.merge(probelinks, udp_sample, how='left')
    probe_to_pr = probe_to_pr.groupby(['link'], sort=False)['pr'].max()
    pr_dict = probe_to_pr.to_dict()
    return pr_dict


def calc_link_to_gene_to_pr(udp, sample_num):
    probelinks = pd.read_csv('data/probelinks.txt', delimiter='\t', header=None)
    probelinks.drop(2, inplace=True, axis=1)
    probelinks.columns = ['probe', 'link']
    probelinks['link'] = probelinks.link.replace('---', 0)
    probelinks['link'] = pd.to_numeric(probelinks.link)
    udp_sample = udp.iloc[:, sample_num].copy().reset_index()
    udp_sample.columns = ['gene', 'pr']
    udp_sample['gene'] = udp_sample.gene.apply(lambda x: x.split('|')[0]) #Remove Entre gene IDs.
    #probelinks['gene'] = probelinks.probe.apply(lambda x: probe_to_gene_dict.get(x, None))
    probelinks = pd.merge(probelinks, probe_to_gene_df, on='probe', how='left')
    probe_to_pr = pd.merge(probelinks, udp_sample, how='left')
    probe_to_pr = probe_to_pr.groupby(['link'], sort=False)['pr'].max()
    pr_dict = probe_to_pr.to_dict()
    return pr_dict


def parse_link(link, link_to_pr_dict, up2ll_dict):
    if ('LL' in link):
        return (link_to_pr_dict.get(int(link.replace('LL:', '')), np.NaN))
    if ('UP' in link):
        return (link_to_pr_dict.get(int(up2ll_dict.get(link.replace('UP:', ''), 0)), np.NaN))
    return np.NaN


#Parse complexes in order to build a dictionary for complex->UDP (similar to the links to UDP disctionary).
#Problem is that complexes might be composed of other complexes. Hence we need another join operation.
#We call complexes that have only proteins (i.e. no other complexes) 'basics'. They are part of other (non-basic) complexes.
#First we calculate UDP for basics, by mutiplying UDPs of corresponding proteins. Then we can calculate
#the UDP of all complexes by, again, multiplying all UDPs of mulecules (proteins and basics) in a path.
#Molecules with no links (i.e. basic complexes) will get a default link. In that case they will get the default UDP value of 1.
def calc_cmplx_to_pr(link_to_pr_dict, up2ll_dict):
    cmplx = pd.read_csv('data/pathologist.complexes.txt', delimiter='\t', header=None)
    cmplx.drop([2, 5, 6, 7], inplace=True, axis=1) # Columns: (0==complex ID, 1==molecule type, 3==molecule ID, 4==Links).

    #Handle molecules (proteins, rnas and compouns). Parse a link string in the links column in the pathways file or the complexes file.
    cmplx['pr'] = cmplx[4].apply(lambda x: max([parse_link(i, link_to_pr_dict, up2ll_dict) for i in str(x).split(',')]))
    cmplx.loc[cmplx[1] == 'rna', 'pr'] = 1
    cmplx.loc[cmplx[1] == 'compound', 'pr'] = 1 # Compounds are assumed to always be present (UDP == 1).

    #Handle basic complexes. The product of molecules on a path are now the accurate basic complex UDP.
    cmplx_basic = cmplx.copy().dropna(subset=[4]) #non-basic complexes don't have links.
    cmplx_basic = cmplx_basic[[0, 'pr']].groupby(0, as_index=False).prod()
    cmplx_basic = cmplx_basic.set_index(0)
    basic_to_dict = cmplx_basic['pr'].to_dict()

    #Handle non-basic complexes (they can be found by searching a molecule of type 'complex').
    #To handle non-basic complex we just need to mutiply the pr for their group of molecules.
    #Update complx table with the basic complexes values, only for molecules of type 'complex' (this will ignore basic-complex since their column 3 has proteins).
    cmplx.loc[cmplx[1] == 'complex', 'pr'] = cmplx[3].apply(lambda x: basic_to_dict.get(x, np.NaN))
    cmplx = cmplx.dropna() #Remove molecules that we couldn't calculate pr for (missing links, unknown basic complexes etc.).
    cmplx = cmplx[[0, 'pr']].groupby(0, as_index=False).prod()
    cmplx = cmplx.set_index(0)
    return cmplx['pr'].to_dict()


def calc_activity_consistency(udp_file, rnaseq_file):
    global probe_to_gene_df
    udp = pd.read_csv(udp_file, index_col=0)
    print(time.ctime(), 'Calculate activity and consistency...')
    orig_paths = pd.read_csv('data/pathologist.db.txt', delimiter='\t', header=None)
    orig_paths.columns = ['path_name', 'path_id', 'molType', 'molName', 'molNum', 'molLink', 'c7', 'c8', 'c9', 'molRole', 'intID', 'intType']
    orig_paths.drop(['molName', 'c7', 'c8', 'c9', 'intType'], inplace=True, axis=1)
    orig_paths = orig_paths.sort_values(['path_id', 'intID'])
    up2ll = pd.read_csv('data/UP2LL.txt', delimiter='\t', header=None, index_col=0)
    up2ll_dict = up2ll[1].to_dict()

    #For RNASeq data we need to create a probe to gene dictionary using Affymetrix table.
    probe_to_gene_df = pd.read_csv('data/HG-U133_Plus_2.na36.annot.csv', index_col=0)
    #Split 'Gene Symbol' in case of multiple genes and stack them in a series.
    new_col = probe_to_gene_df['GeneSymbol'].str.split('|', expand=True).stack()
    #Join the new series (after setting a name and index) back to the dataframe.
    probe_to_gene_df = probe_to_gene_df.join(pd.Series(index=new_col.index.droplevel(1), data=new_col.values, name = 'gene'))
    probe_to_gene_df.drop_duplicates(inplace=True)
    del probe_to_gene_df['GeneSymbol']
    probe_to_gene_df.reset_index(inplace=True)
    probe_to_gene_df.columns = ['probe', 'gene']

    sample_results = np.empty((0,6)) # Table columns are (sample_num, path_id, activity, consistency)
    for sample_num in range(0, len(udp.columns)):
        if (sample_num % 10 == 0):
           print(f'Processing sample {sample_num}...')
           #sys.stdout.flush()
        #Calculate UDP of the molecules in all the paths.
        if(rnaseq_file == 0):
            link_to_pr_dict = calc_link_to_pr(udp, sample_num)
        else:
            link_to_pr_dict = calc_link_to_gene_to_pr(udp, sample_num)
        cmplx_to_pr_dict = calc_cmplx_to_pr(link_to_pr_dict, up2ll_dict)
        paths = orig_paths.copy()
        #gc.collect()
        #If a molecule does not have a probability we need to remove the whole interaction from activity and consistency calculations.
        paths.loc[paths.molType == 'protein', 'pr'] = paths.molLink.apply(lambda x: max([parse_link(i, link_to_pr_dict, up2ll_dict) for i in str(x).split(',')])) #Max returns NaN if there is at least one NaN in the list.
        paths.loc[paths.molType == 'compound', 'pr'] = 1 #Compounds are assumed to always be present
        paths.loc[paths.molType == 'rna', 'pr'] = 1
        paths.loc[paths.molType == 'complex', 'pr'] = paths.molNum.apply(lambda x: cmplx_to_pr_dict.get(x, np.NaN))

        #Calculate activity and consistency of each interaction.
        #Activity of interactions is a multiplication of the inputs.
        paths.loc[paths.molRole == 'inhibitor', 'pr'] = 1 - paths.pr
        # max(paths['path_id'].groupby(paths['intID']).nunique()) #This shows that some interactions participate in few paths.
        # Save the pr for output molecules, since it also receieved a UDP.
        paths['pr_output'] = 1
        paths.loc[paths.molRole == 'output', 'pr_output'] = paths.pr
        paths.loc[paths.molRole == 'output', 'pr'] = 1
        #Handle molecules that we could not find UDP for, we need to remove the whole interaction.
        #Pandas like R is ignoring NA values in gropuby.transform, so instead we zero it.
        #paths = paths.dropna(subset=['pr'])
        paths = paths.fillna({'pr':0})
        paths['interaction_activity'] = paths.groupby(['path_id', 'intID'])['pr'].transform('prod')
        paths['interaction_consistency'] = paths.interaction_activity * paths.pr_output + (1 - paths.interaction_activity) * (1 - paths.pr_output)

        #Calculate activity and consistency of each path.
        #Both activity and consistency are averages of all interactions activities and consistencies.
        paths['activity'] = paths.groupby('path_id')['interaction_activity'].transform('mean')
        paths['consistency'] = paths.groupby('path_id')['interaction_consistency'].transform('mean')
        paths['sample_num'] = sample_num
        paths_result = paths[['sample_num', 'path_name', 'path_id', 'activity', 'consistency', 'molRole']].drop_duplicates()
        sample_results = np.concatenate((sample_results, paths_result.values)) # axis=0
    results_df = pd.DataFrame(data=sample_results, columns=['sampleID', 'path_name', 'pathID', 'Activity', 'Consistency', 'molRole'])
    sample_names_dict = dict(enumerate(udp.columns))
    results_df['sample_name'] = results_df.sampleID.apply(lambda x: sample_names_dict.get(x, 'Unknown'))
    results_df.to_csv('data/output_activity.csv', index=False)
    print(time.ctime(),"Done.")
    return results_df, paths