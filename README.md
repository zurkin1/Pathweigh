# PathWeigh A Biochemical Pathway Analysis Tool
![Pathway analysis](https://norbis.w.uib.no/files/2016/05/F1.large_-768x623.jpg)

## Overview
PathWeigh is a Python-based pathway analysis tool tailored for microarray and RNA-seq data analysis. It employs a unique graph-based algorithm to enable the analysis of diverse cellular states, such as T cell subtypes. Designed to be open-source, extensible, and computationally efficient, PathWeigh provides researchers with a versatile framework for uncovering biologically meaningful insights from high-dimensional transcriptomics data.

## Key Features
- Tailored for microarray and RNA-seq data analysis
- Unique graph-based algorithm for pathway analysis
- Efficient classification of diverse cellular states
- Open-source and extensible
- Computationally efficient

## Installation
Simply clone this repository using git clone command.

## Quick Start
```
df = pd.read_csv('./data/input.csv')
calc_activity(df)
```

## Usage
For detailed usage instructions, please refer to the notebooks in the `notebooks` folder, which demonstrate RNA-seq data processing using PathWeigh.

## Supported Pathways
PathWeigh currently supports 581 curated pathways. Click the link to view the full list. [List of supported pathways.](data/pathnames.txt)
[Guide for adding a new pathway.](data/guide.md)

## Contributing
We welcome contributions! Please see our Contributing Guidelines for more information on how to get involved.

## License
PathWeigh is available under the MIT license. See the LICENSE file for more details.

## Support
For questions, issues, or feature requests, please open an issue on our GitHub repository.

For additional support, contact: zurkin at yahoo dot com.

## Citation
If you use PathSingle in your research, please cite our paper:

Livne, D., Efroni, S. (2022). PathWeigh – Quantifying the Behavior of Biochemical Pathway Cascades. In: Rojas, I., Valenzuela, O., Rojas, F., Herrera, L.J., Ortuño, F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2022. Lecture Notes in Computer Science(), vol 13347. Springer, Cham. https://doi.org/10.1007/978-3-031-07802-6_29

## Acknowledgments
We thank the scientific community for their valuable feedback and contributions to this project.