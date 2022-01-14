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

To begin the tutorial, please unzip the tutorial's expression data.

```bash
$ cd CellLayers
$ unzip /Data/*/\*.zip 
```
---

To recreate **Figure 1A**, please follow [PBMC_Tutorial.ipynb](https://github.com/apblair/CellLayers/blob/master/Notebooks/PBMC_Tutorial.ipynb).

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
<!-- <image src="Images/example.png"> -->
    
**Fig. 1 (A)** PBMC multi-resolution analysis from 0.1 to 0.5. Edges are painted by _CD3E_, which is a marker gene for CD8+ T, Memory CD4 T, and Naive CD4 T cells. Nodes are painted by Silhouette score. The lower Sihouette values indicate samples are near the decision boundary of neighboring clusters.

---
    
To recreate **Figure 1B**, please follow [PBMC_Tutorial.ipynb](https://github.com/apblair/CellLayers/blob/master/Notebooks/PBMC_Tutorial.ipynb).

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
    
**Fig. 1 (B)** Nodes painted by enrichR GO 2018 Biological Process gene set scores for GO:0002480. The node hovertemplate provides users cluster performance metrics (modularity and silhouette scores), GO term, enrichR gene set score, and the top 5 differentially expressed genes. Edges are colored by NK marker gene _CD8A_.
    
---
    
To recreate **Figure 1C**, please follow [TBX5_Tutorial.ipynb](https://github.com/apblair/CellLayers/blob/master/Notebooks/TBX5_Tutorial.ipynb):
    
**Fig. 1 (C)** iPSC-derived cardiomyocyte multi resolution analysisfrom 0.1 to 0.5. Edges are painted by coexpression of _TNNT2_ (red), _COL1A1_ (green), and _NR2F2_ (blue). Nodes are painted by Silhouette score. Arrows on the Ternary plot indicate the direction of the co-expression scale for each edge in the Sankey chart.
  
## Tutorial: Generate input data using SetupCellLayers

To create the PBMC CellLayers input data using **SetupCellLayers**, please follow [PBMC_SetupCellLayers.ipynb](https://github.com/apblair/CellLayers/blob/master/Notebooks/PBMC_SetupCellLayers.ipynb).

```R
library(enrichR)
library(Seurat)
library(SetupCellLayers)

sobj <- readRDS('CellLayers/Data/PBMC/pbmc3k_CellLayers.rds')

cl.setup <- compute_modularity(sobj, seq(0.1, 0.5, by=0.1),
                               'CellLayers/Data/Test-SetupCellLayers/PBMC_modularity.csv')

cl.setup[['sil']] <- compute_silhouette_scores(cl.setup$sobj,
                                               as.character(cl.setup$mod$'resolution'),
                                               'CellLayers/Data/Test-SetupCellLayers/PBMC_silhouette_scores.csv')

compute_enrichment("GO_Biological_Process_2018", 
                   cl.setup$sobj, 
                   seq(0.1, 0.5, by=0.1), 
                   'CellLayers/Data/Test-SetupCellLayers/')

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


