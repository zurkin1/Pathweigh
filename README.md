# PathWeigh - A Biochemical Pathway Analysis Tool
![Pathway analysis](https://norbis.w.uib.no/files/2016/05/F1.large_-768x623.jpg)

## The Concept of Biological Pathway
A biological pathway is a series of molecular interactions and reactions that govern a particular cellular process or function [74]. It represents a coordinated sequence of biochemical events involving molecules such as proteins, genes, metabolites, and other biomolecules. These molecular components interact with each other in a specific order and regulate each other's activities, ultimately leading to a specific biological outcome or phenotype.

Biological pathways can be thought of as complex networks or circuits, where each component (gene, protein, metabolite) acts as a node, and the interactions or reactions between them form the connections or edges. These connections can represent various types of relationships, such as enzyme-substrate interactions, gene regulation, protein modifications, or signal transduction events.

Pathways are typically organized into distinct categories based on the cellular processes they govern, such as: metabolic pathways, signaling pathways, gene regulation pathways, disease pathways.

Biological pathways are typically represented using graphical notations or computational models, where nodes represent the molecular components, and edges represent the interactions or reactions between them. These representations capture the intricate relationships and regulatory mechanisms that govern cellular functions and allow researchers to study and understand the complex dynamics of biological systems.

## Pathway Analysis
Pathway analysis is a powerful approach in systems biology and bioinformatics that aims to understand biological processes and disease mechanisms by studying the coordinated interactions and relationships among genes, proteins, and metabolites within the context of molecular pathways. Rather than focusing on individual molecules or genes in isolation, pathway analysis considers the collective behavior of interconnected components that work together to carry out specific biological functions. At its core, pathway analysis leverages prior knowledge of biological pathways, which are networks that represent the complex interactions and regulatory mechanisms underlying various cellular processes. These pathways are typically curated from existing literature, databases, and experimental evidence, and they provide a structured framework for interpreting and contextualizing large-scale molecular data, such as gene expression profiles, protein abundances, or metabolite levels.

The primary goal of pathway analysis is to identify pathways or sets of pathways that are significantly dysregulated or perturbed given a biological condition or disease state, compared to a control or reference state. By connecting molecular data to known biological pathways, researchers can gain insights into the functional implications of observed changes and identify potential mechanisms driving disease progression or response to therapeutic interventions.

Some methods determine whether a set of differentially expressed genes, proteins, or metabolites is enriched or over-represented in specific pathways compared to what would be expected by chance.

Other methods use pathway-level statistics to compute and aggregate the evidence from individual genes or molecules within a pathway. This statistic can then be used to rank pathways based on their overall dysregulation or perturbation. That is the approach taken in our work.

Topology-based methods: These methods consider the topology or structure of the pathways, accounting for the interactions and dependencies among pathway components. Examples include gene set enrichment analysis (GSEA) and impact analysis, which incorporate information about the position and relationships of genes or molecules within the pathway.

Pathway analysis has several advantages over traditional gene-centric or molecule-centric approaches. First, it provides a more biologically interpretable and meaningful context for understanding molecular data, as pathways represent well-characterized biological processes. Second, by considering the collective behavior of multiple genes or molecules, pathway analysis can reveal subtle but coordinated changes that may be missed when analyzing individual components separately. Third, pathway analysis can help identify potential therapeutic targets or biomarkers by highlighting dysregulated pathways that are critical for disease progression or response to treatment.

## Running PathWeigh
In order to run PathWeigh please refer to the notebook folder where twe notebooks are availabe. One demonstrating bulk RNAseq data, and the other single-cell RNAseq data.


### [List of supported pathways.](data/pathnames.txt)

### [Guide for adding a new pathway.](data/guide.md)

##
Livne, D., Efroni, S. (2022). PathWeigh – Quantifying the Behavior of Biochemical Pathway Cascades. In: Rojas, I., Valenzuela, O., Rojas, F., Herrera, L.J., Ortuño, F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2022. Lecture Notes in Computer Science(), vol 13347. Springer, Cham. https://doi.org/10.1007/978-3-031-07802-6_29.