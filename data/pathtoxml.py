import numpy as np
import pandas as pd
from activity import path_activity

colors_green_to_red = ['00FF00','11FF00','22FF00','33FF00','44FF00','55FF00','66FF00','77FF00','88FF00',
          '99FF00','AAFF00','BBFF00','CCFF00','DDFF00','EEFF00','FFFF00','FFEE00','FFDD00',
          'FFCC00','FFBB00','FFAA00','FF9900','FF8800','FF7700','FF6600','FF5500','FF4400',
          'FF3300','FF2200','FF1100','FF0000']


def print_and_wait(some_long_message):
    lines = some_long_message.split('\n')
    i=0
    while i < len(lines):
        print ('\n'.join(lines[i:i+10]))
        inp = input("press enter to read more, q to quit...")
        if(inp == 'q'):
            return
        i += 10


"""
    <entry id="29" name="hsa:4967 hsa:{55753}" type="protein" 
        link="https://www.kegg.jp/dbget-bin/www_bget?hsa:4967+hsa:{55753}">
        <graphics name="{OGDH(0.9)}" fgcolor="#000000" bgcolor="#BFFFBF"
             type="rectangle" x="661" y="579" width="46" height="17"/>
    </entry>
"""
def xmlparser(path_id, sample_num, udp_file='/tmp/output_udp.csv', rnaseq=0):
    nodes = set()
    intxns = set()

    paths = pd.read_csv('pathologist.db.txt', delimiter='\t', header=None)
    paths.columns = ['path_name', 'path_id', 'molType', 'molName', 'molNum', 'molLink', 'c7', 'c8', 'c9', 'molRole', 'intID', 'intType']
    print(f'Create Kegg XML for path: {path_id}, sample: {sample_num}')
    if (isinstance(path_id, int)):
        paths = paths.loc[paths.path_id == path_id].copy()
    else:
        #If instead of path_id the input has the path name.
        paths = paths.loc[paths.path_name == path_id]
        path_id = paths.iloc[0]['path_id']
    #Add '+' in case of an active molecule.
    paths.loc[paths.c7 == 'active', 'molName'] = paths.molName + '+'

    #We need to get the extra information of interactions activity (that is not saved in the standard calc_activity_consistency run).
    activity_obj = path_activity(udp_file, rnaseq)
    sampleID = activity_obj.udp.columns[sample_num]
    _, prbs = activity_obj.process_samples([sampleID])
    prbs = prbs.loc[prbs['path_id'] == path_id]
    prbs = prbs[['molNum', 'molRole', 'pr', 'pr_output', 'intID', 'interaction_activity', 'interaction_consistency']]
    paths = pd.merge(paths, prbs, on=['molNum', 'intID', 'molRole'], how='left')

    xml = '<?xml version="1.0"?>' +\
          '<!DOCTYPE pathway SYSTEM "https://www.kegg.jp/kegg/xml/KGML_v0.7.2_.dtd">' +\
          '<pathway name="path:hsa{0}" org="hsa" number="{1}" title="{2}">'.format(paths.iloc[0,1], paths.iloc[0,0], paths.iloc[0,1])

    for index, row in paths.iterrows():
        if(isinstance(row['molName'], float)): #Handle cases of empty molecule name.
            row['molName'] = str(row['molNum'])
        if (row['molName'] not in nodes):
            nodes.add(row['molName'])

            bg_color = "42ECEF" #Light turquoise
            if (row['molRole'] == 'output'):
                fg_color = int(np.nan_to_num(row["pr_output"])*30)
                pr = str(np.nan_to_num(row["pr_output"]))
            else:
                fg_color = int(np.nan_to_num(row["pr"])*30)
                pr = str(np.nan_to_num(row["pr"]))

            if (row['molType'] == 'complex'):
                row['molType'] = 'ortholog' #Kegg uses 'group' for complexes.
                bg_color = "C3B6B6" #Grey
            if (row['molType'] == 'compound'):
                row['molType'] = 'protein'

            xml += f'<entry id="{row["molName"]}" name="hsa:4967 hsa:{row["molNum"]}" type="{row["molType"]}" ' +\
                   f'link="https://www.kegg.jp/dbget-bin/www_bget?hsa:4967+hsa:{row["molNum"]}"> ' +\
                   f'<graphics name="{row["molName"]}" fgcolor="#{colors_green_to_red[fg_color]}" bgcolor="#{bg_color}" x="{pr}" ' +\
                   f'type="rectangle" width="{min(60, 8*len(row["molName"]))}" height="17"/> ' +\
                   '</entry>'
        #Add the interaction
        if(row['intID'] not in intxns):
            intxns.add(row['intID'])
            bg_color = "FFF8F8" #White
            fg_color = int(row["interaction_activity"]*30)
            pr_act = str(row["interaction_activity"])
            pr_con = str(row["interaction_consistency"])
            # Add the junction representing the interaction.
            xml += f'<entry id="{row["intID"]}" name="{row["intID"]}" type="reaction" ' +\
                   'link="https://www.kegg.jp/dbget-bin/www_bget?C00158"> ' +\
                   f'<graphics name="{row["intType"]}" fgcolor="#{colors_green_to_red[fg_color]}" bgcolor="#{bg_color}" x="{pr_act}" y="{pr_con}" ' +\
                   f'type="circle" width="{min(60, 8*len(row["intType"]))}" height="15"/> ' +\
                   '</entry>'

        #Add the graph edges.
        if(row['molRole'] == 'input' or row['molRole'] == 'agent'):
            xml += f'<relation entry1="{row["molName"]}" entry2="{row["intID"]}" type="PPrel"> '\
                        '<subtype name="activation" value="-->"/> '\
                   '</relation>'
        elif(row['molRole'] == 'output'):
            xml += f'<relation entry1="{row["intID"]}" entry2="{row["molName"]}" type="PPrel"> '\
                        '<subtype name="activation" value="-->"/> '\
                   '</relation>'
        elif(row['molRole'] == 'inhibitor'):
            xml += f'<relation entry1="{row["molName"]}" entry2="{row["intID"]}" type="PPrel"> '\
                        '<subtype name="activation" value="--|"/> '\
                   '</relation>'
    xml += '</pathway>'
    text_file = open(f"data/output_path_{path_id}_{sample_num}.xml", "w")
    text_file.write(xml)
    text_file.close()
    print(f'Output XML file: data/output_path_{path_id}_{sample_num}.xml')


if __name__ == "__main__":
    paths = pd.read_csv('data/pathologist.db.txt', delimiter='\t', header=None)
    paths.columns = ['path_name', 'path_id', 'molType', 'molName', 'molNum', 'molLink', 'c7', 'c8', 'c9', 'molRole', 'intID', 'intType']
    path_id = 0 #100151 #726290
    sample_id = 0
    print('########## Pathway to Kegg KGML ##########')
    while (path_id == 0 or path_id =='l'):
        path_id = input('Select a path ID (l to list paths):')
        if (path_id == 'l'):
            print_and_wait(paths[['path_name','path_id']].drop_duplicates().to_string(index = False))
    path_id = int(path_id)

    samples = pd.read_csv('data/output_activity.csv')
    sample_id = -1
    while (sample_id == -1 or sample_id =='l'):
        sample_id = input('Select a sample ID (l to list samples):')
        if(sample_id == 'l'):
            print_and_wait(samples[['sampleID', 'sample_name']].drop_duplicates().to_string(index = False))
    sample_id = int(sample_id)

    xmlparser(paths, path_id, sample_id, 'data/output_udp_acc.csv', 1)
    print ('########## Done.')