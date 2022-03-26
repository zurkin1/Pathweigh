import numpy as np
import pandas as pd
import time
import pickle


def process_samples(udp_chunk):
    sample_results = np.empty((0, 6))
    mapping = kegg['mapping']
    #For every sample Calculate activity and consistency of a path. Output table columns are (sampleID, path_id, activity, consistency)
    for sample in udp_chunk:
      path_inter_ids = kegg['pathList'][path_id]
      for inter_id in path_inter_ids:
        interaction = kegg['interactionList'].loc[kegg['interactionList']['interaction.id'] == inter_id]
        inputm = interaction.iloc[0]['input']
        outputm = interaction.iloc[0]['output']
        node_in = mapping.loc[mapping['node.id'] == inputm]
        node_out = mapping.loc[mapping['node.id'] == outputm]
        try:
          activity_in = udp[sample].loc[node_in['symbol'].values].prod()
          activity_out = udp[sample].loc[node_out['symbol'].values].prod()
        except:
          print('.', end="")

        #nn = kegg['node.name']
        #nt = kegg['node.type']
        #ind = nn.index('RPS6KB2')
        #nt[ind]

'''
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
      # paths = paths.dropna(subset=['pr'])
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
'''


udp = pd.read_csv('Pathweigh/data/output_udp_zurkin.csv', index_col = 0)
with open('pid', 'rb') as f:
  paths = pickle.load(f)
kegg = paths['PID.db']['KEGG']
path_id = 'hsa05221'
process_samples(udp)
# activity_obj.calc_activity_consistency_multi_process()
# activity_obj.graphparser(100035, 1)