import numpy as np
import pandas as pd
import time
import pickle


def process_samples(udp_chunk):
    mapping = kegg['mapping']
    path_inter_ids = kegg['pathList'][path_id]
    for inter_id in path_inter_ids:
        interaction = kegg['interactionList'].loc[kegg['interactionList']['interaction.id'] == inter_id]
        inputm = interaction.iloc[0]['input']
        outputm = interaction.iloc[0]['output']
        node_in = mapping.loc[mapping['node.id'] == inputm]
        node_out = mapping.loc[mapping['node.id'] == outputm]
        total_activity = 0
        total_inter = 0
        try:
          activity_in = udp[sample].loc[[str.lower(x) for x in node_in['symbol'].values]].prod()
          activity_out = udp[sample].loc[[str.lower(x) for x in node_out['symbol'].values]].prod()
          total_activity += activity_in
          total_inter += 1
        except:
          print('.', end="")

    #nn = kegg['node.name']
    #nt = kegg['node.type']
    #ind = nn.index('RPS6KB2')
    #nt[ind]


#udp = pd.read_csv('Pathweigh/data/output_udp_zurkin.csv', index_col = 0)
with open('./data/pid', 'rb') as f:
  paths = pickle.load(f)
for i in paths['PID.db']['KEGG']:
  print(i, type(paths['PID.db']['KEGG'][i]), len(paths['PID.db']['KEGG'][i]))
kegg = paths['PID.db']['KEGG']
#path_id = 'hsa05221'
#process_samples(udp)
pl = kegg['pathList']
il = kegg['interactionList']
mapping = kegg['mapping']
nn = kegg['node.name']
nt = kegg['node.type']

#pd.DataFrame.from_dict(pl, orient='index').T.to_csv('pathlist.csv')

#k = pd.DataFrame({'node_name': nn, 'node_type': nt})
#k.to_csv('node.csv')

#il.to_csv('interactions.csv')

mapping.to_csv('mappings.csv')