# Pathweigh A Biochemical Pathway Analysis Tool
![Pathway analysis](https://norbis.w.uib.no/files/2016/05/F1.large_-768x623.jpg)

## Requirements

- numpy
scipy
pandas
scikit-learn
ipywidgets
networkx
plotly
bs4.

## Running
- Input: a csv file called input.csv that contains samples and gene RNA expression levels.
- All functions are available in udp.py.
- Python code:
```
    #Calculate UDP.
    #data = pd.read_csv(relative_path + './input.csv', index_col=0)
    #Comparison of various fitting methods.
    for func in [calc_udp_poisson, calc_udp_nbm, calc_udp_gmm, calc_udp_norm, calc_udp_gennorm]:
        udp, aic = func(sample_data, aic_test=True)
        print(f'Function: {func.__name__}, aic: {aic}')
    calc_udp_multi_process(data, True)

    #Calculate activity.
    udp = pd.read_csv('./data/output_udp.csv', index_col = 0)
    with open('./data/pid', 'rb') as f:
        paths = pickle.load(f)
    kegg = paths['PID.db']['KEGG']
    path_id = 'hsa05221'
    sample = udp['17-002']
    process_sample(kegg, path_id, sample)
```
## Pathway IDs
https://github.com/NCIP/pathway-interaction-database/tree/master/download.
