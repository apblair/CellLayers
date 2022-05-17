![BuildStatus](https://github.com/apblair/CellLayers/actions/workflows/main.yml/badge.svg?event=push)

# Cell Layers: Uncovering clustering structure in unsupervised single-cell transcriptomic analysis

Cell Layers is an interactive Sankey tool for the quantitative investigation of gene expression, coexpression, biological processes, and cluster integrity across clustering resolutions.

## Installation
The **CellLayers** Python module can be installed via pip after cloning the repository. Please note, the Python module requires the cluster metric and enrichment data to be generated independently or via our R library **SetupCellLayers**. 

```bash
$ git clone git@github.com:apblair/CellLayers.git
$ cd CellLayers
$ pip install .
```

The **SetupCellLayers** R library can be installed via devtools after cloning the repository.
```R
> library(devtools)
> setwd("CellLayers")
> devtools::install_github("apblair/CellLayers/SetupCellLayers")
> library(SetupCellLayers)
```
## Documentation
Please consider citing Cell Layers if you used the application or code snippets for your analysis.
    
**Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis** Andrew P. Blair, Robert K. Hu, Elie N. Farah, Neil C. Chi, Katherine S. Pollard, Pawel F. Przytycki*, Irfan S. Kathiriya*, Benoit G. Bruneau*
bioRxiv 2020.11.29.400614; doi: https://doi.org/10.1101/2020.11.29.400614
    
    
## Authors Contributions
A.P.B. conceived and initiated the project. R.H assisted in analysis. B.G.B., P.F.P, and I.S.K. supervised A.P.B. I.S.K. provided datasets. K.S.P. advised. All authors commented on the manuscript.

## Acknowledgments
We thank the Cytoscape and scNetViz developers Alex Pico and Scooter Morris for their input on Plotly, Dan Carlin for his input on multi-resolution analysis, members of the CIRM Heart of Cells group, Gladstone Bioinformatics core, and Bruneau lab for discussions and comments. 

## Funding
California Institute for Regenerative Medicine (RB4-05901 to B.G.B)


