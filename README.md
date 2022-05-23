![BuildStatus](https://github.com/apblair/CellLayers/actions/workflows/main.yml/badge.svg?event=push)

# Cell Layers: Uncovering clustering structure in unsupervised single-cell transcriptomic analysis

Cell Layers is an interactive Sankey tool for the quantitative investigation of gene expression, coexpression, biological processes, and cluster integrity across clustering resolutions. 

For more information, please visit our [Read the Docs](https://celllayers.readthedocs.io/en/latest/index.html). 

![plot](/images/example.png)

**Keywords**: single cell, clustering, gene expression, visualization, data analysis

## Installation
The **CellLayers** Python module can be installed via pip. Please note, the Python module requires the cluster metric and enrichment data to be generated independently or via our R library **SetupCellLayers**. 

```bash
$ pip install CellLayers
```

The **SetupCellLayers** R library can be installed via devtools.
```R
> library(devtools)
> options(timeout=9999999)
> devtools::install_github("apblair/CellLayers/SetupCellLayers")
> library(SetupCellLayers)
```

Both **CellLayers** and **SetupCellLayers** are fully containerized via Docker and Singularity, which are extended from the jupyter/datascience-notebook image.
```bash
$ docker pull apblair/cell-layers:v0.2
$ docker run -it --rm -p 10000:8888 -v "${PWD}":/home/jovyan/work apblair/cell-layers:v0.2
```

```bash
$ singularity pull --arch amd64 library://apblair/single-cell-tools/cell-layers:v0-2
$ singularity exec cell-layers_v0-2.sif start.sh jupyter lab --port=9595
```

## Tutorial

Please see the [notebooks](https://github.com/apblair/CellLayers/tree/master/notebooks) folder for a tutorial on how to reproduce Fig.1 using **SetupCellLayers** and **CellLayers**.

## Reference
    
**Cell Layers: Uncovering clustering structure in unsupervised single-cell transcriptomic analysis** Andrew P. Blair, Robert K. Hu, Elie N. Farah, Neil C. Chi, Katherine S. Pollard, Pawel F. Przytycki*, Irfan S. Kathiriya*, Benoit G. Bruneau*
bioRxiv 2020.11.29.400614; doi: https://doi.org/10.1101/2020.11.29.400614

Please consider citing Cell Layers if you used the application or it's code snippets for your project.
    
## Authors Contributions
A.P.B. conceived and initiated the project. R.H assisted in analysis. B.G.B., P.F.P, and I.S.K. supervised A.P.B. I.S.K. provided datasets. K.S.P. advised. All authors commented on the manuscript.

## Acknowledgments
We thank the Cytoscape and scNetViz developers Alex Pico and Scooter Morris for their input on Plotly, Dan Carlin for his input on multi-resolution analysis, members of the CIRM Heart of Cells group, Gladstone Bioinformatics core, and Bruneau lab for discussions and comments. 

## Funding
California Institute for Regenerative Medicine (RB4-05901 to B.G.B), Additional Ventures (to B.G.B and K.S.P), National Institutes of Health (5U01HL157989-02 to B.G.B and K.S.P), UCSF Department of Anesthesia and Perioperative Care (to I.S.K.), and the Younger family fund (to B.G.B).

