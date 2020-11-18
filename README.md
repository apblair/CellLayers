# Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis
#### Andrew P. Blair, Robert K. Hu, Katherine S. Pollard, Pawel F. Przytycki*, Irfan S. Kathiriya*, Benoit G. Bruneau*
*Contact*: andrew.blair@gladstone.ucsf.edu and irfan.kathiriya@ucsf.edu
<image src="Examples/example.png">

## Motivation
Unsupervised clustering of single-cell transcriptomics is a powerful method for identifying cell populations. Static visualization techniques for single-cell clustering only display results for a single resolution parameter. Analysts will often evaluate more than one resolution parameter, but then only report one.

## Results
We developed Cell Layers, an interactive Sankey tool for the quantitative investigation of gene expression, coexpression, biological processes, and cluster integrity across clustering resolutions. Cell Layers enhances the interpretability of single-cell clustering by linking molecular data and cluster evaluation metrics to provide novel insight into cell populations.

## Installation
Cell Layers may be installed via pip, conda, or Docker.

```bash
$ pip install CellLayers
```
```bash
$ conda install CellLayers
```
```bash
$ docker pull CellLayers
```

## Running Container

```bash
$ docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work apblair/computing-envs:ab-JupyterLab_seurat_v3.2.0_scanpy_v1.5.0
```

## Seurat Data Generation
```R
library(Seurat)
```

## Usage
```python
import CellLayers
import pandas as pd

pbmc_exp = pd.read_csv('Data/PBMC_exp.csv', index_col=[0]) # cell by gene expression matrix
pbmc_meta = pd.read_csv('Data/PBMC_meta.csv', index_col=[0]) # cell by resolution matrix
pbmc_modularity = pd.read_csv('Data/cluster_metrics/pbmc_modularity.csv', index_col=[0])
pbmc_sil = pd.read_csv('Data/cluster_metrics/pbmc_silhouette_scores.csv', index_col=[0])

sankey_dict = CellLayers(pbmc_exp, 
                         pbmc_meta, 
                         modularity=pbmc_modularity,
                         silhouette=pbmc_sil,
                         genes=['MS4A1', 'GNLY',
                                'CD3E', 'FCER1A'],
                         tri_coexpressed_genes=[['MS4A1',
                                                 'FCER1A', 
                                                 'FCGR3A']])
```



## Documentation
Please consider citing Cell Layers if you used the application for your analysis.

## Authors Contributions
A.P.B. conceived and initiated the project. R.H assisted in analysis. B.G.B., P.F.P, and I.S.K. supervised A.P.B. I.S.K. provided datasets. K.S.P. advised. All authors commented on the manuscript.

## Acknowledgments
We thank the Cytoscape and scNetViz developers Alex Pico and Scooter Morris for their input on Plotly, Dan Carlin for his input on multi-resolution analysis, members of the CIRM Heart of Cells group, Gladstone Bioinformatics core, and Bruneau lab for discussions and comments. 

## Funding
California Institute for Regenerative Medicine (RB4-05901 to B.G.B)


