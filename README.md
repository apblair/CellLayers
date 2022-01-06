# Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis

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

We have also provided a Docker image, which encapsulates the environment to run **CellLayers** and **SetupCellLayers**. 
```bash
$ docker pull bruneaulab/cell-layers:0.1
```

## Running the Container

```bash
$ docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work bruneaulab/cell-layers:0.1
```

For HPC clusters that support Linux containers, the Docker container can be pulled and built as a Singularity container. After building the Singularity container, an interactive session can be run with these commands:
```bash
$ singularity build cell-layers.sif docker://bruneaulab/cell-layers:0.1
$ singularity exec cell-layers.sif start.sh jupyter lab --port=9595
```

## Tutorial: Reproduce publication figure

Please unzip the tutorial's expression data.

```bash
$ cd CellLayers/Data
$ unzip PBMC_exp.csv.zip 
```
---

To recreate **Figure 1A**, please open CellLayers/Notebooks/PBCM_Tutorial.ipynb in a Jupyter environment and run the following cell:

```python
import CellLayers
import pandas as pd

pbmc_exp = pd.read_csv('CellLayers/Data/PBMC_exp.csv', index_col=[0]) # cell by gene expression matrix
pbmc_meta = pd.read_csv('CellLayers/Data/PBMC_meta.csv', index_col=[0]) # cell by resolution matrix
pbmc_modularity = pd.read_csv('CellLayers/Data/pbmc_modularity.csv', index_col=[0]) # cluster resolution modularity scores
pbmc_silhouette_scores = pd.read_csv('CellLayers/Data/pbmc_silhouette_scores.csv', index_col=[0]) # cluster resolution community silhouette scores

sankey_fig, sankey_dict = CellLayers.build_sankey(pbmc_exp,
                                 pbmc_meta,
                                 modularity=pbmc_modularity,
                                 silhouette=pbmc_silhouette_scores,
                                 genes=['CD3E','FCER1A'])
sankey_fig.show()
```
<image src="Images/example.png">
    
**Fig. 1 (A)** PBMC mulit-resolution analysis from 0.1 to 0.5. Edges are painted by CD3E, which is a marker gene for CD8+ T, Memory CD4 T, and Naive CD4 T cells. Nodes are painted by Silhouette score. The lower Sihouette values indicate samples are near the decision boundary of neighboring clusters.

---
    
To recreate **Figure 1B**, please open CellLayers/Notebooks/PBMC_Tutorial.ipynb in a Jupyter environment and run the following cell:

```Python
pbmc_enrichment = pd.read_csv('CellLayers/Data/PBMC/pbmc_enrichment.csv', index_col=[0]) # geneset, cluster resolution communities, and combined score
pbmc_top_genes = pd.read_csv('CellLayers/Data/PBMC/pbmc_top_genes.csv', index_col=[0]) # cluster resolution communities and top genes

# Create a list of your gene set(s) of interest
geneset_oi = ['antigen processing and presentation of exogenous peptide antigen via MHC class II (GO:0019886)',
              'antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-independent (GO:0002480)']

genes = ['CD3E']
enrichment_sankey_fig, enrichment_sankey_dict = CellLayers.build_enrichment_sankey(sankey_dict,
                                                                            geneset_oi,
                                                                            pbmc_enrichment,
                                                                            pbmc_top_genes,
                                                                            genes)
enrichment_sankey_fig.show()
```
    
**Fig. 1 B** Nodes painted byenrichR GO 2018 Biological Process gene set scores for GO:0002480. The node hovertemplate provides users cluster performance metrics (modularity and silhouette scores), GOterm title, enrichR gene set score, and the top 5 differentially expressed genes. Edges arecolored by NK marker gene ​CD8A​.
    
---
    
To recreate **Figure 1C**, please open CellLayers/Notebooks/TBX5_Tutorial.ipynb in a Jupyenfarahter environment and run the following cell:

**Fig. 1 C** iPSC-derived cardiomyocyte multi resolution analysisfrom 0.1 to 0.5. Edges are painted by coexpression of ​TNNT2 ​(red), ​COL1A1 ​(green), andNR2F2 ​(blue). Nodes are painted by Silhouette score. Arrows on the Ternary plot indicate thedirection of the co-expression scale for each edge in the Sankey chart.
  
## Tutorial: Generate input data using SetupCellLayers
    
```R
library(enrichR)
library(Seurat)
library(SetupCellLayers)
sobj <- readRDS('../Data/PBMC/pbmc3k_CellLayers.rds')


```
    
## Documentation
Please consider citing Cell Layers if you used the application or code snippets for your analysis.
    
**Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis** Andrew P. Blair, Robert K. Hu, Elie N. Farah, Neil C. Chi, Katherine S. Pollard, Pawel F. Przytycki, Irfan S. Kathiriya, Benoit G. Bruneau
bioRxiv 2020.11.29.400614; doi: https://doi.org/10.1101/2020.11.29.400614
    
    
## Authors Contributions
A.P.B. conceived and initiated the project. R.H assisted in analysis. B.G.B., P.F.P, and I.S.K. supervised A.P.B. I.S.K. provided datasets. K.S.P. advised. All authors commented on the manuscript.

## Acknowledgments
We thank the Cytoscape and scNetViz developers Alex Pico and Scooter Morris for their input on Plotly, Dan Carlin for his input on multi-resolution analysis, members of the CIRM Heart of Cells group, Gladstone Bioinformatics core, and Bruneau lab for discussions and comments. 

## Funding
California Institute for Regenerative Medicine (RB4-05901 to B.G.B)


