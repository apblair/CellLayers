# Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis
#### Andrew P. Blair, Robert K. Hu, Irfan S. Kathiriya, Katherine S. Pollard, Pawel F. Przytycki, Benoit G. Bruneau

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

## Usage
```python
import CellLayers
import pandas as pd

pbmc_exp = pd.read_csv('Data/PBMC_exp.csv', index_col=[0]) # cell by gene expression matrix
pbmc_meta = pd.read_csv('Data/PBMC_meta.csv', index_col=[0]) # cell by resolution matrix
pbmc_modularity = pd.read_csv('Data/cluster_metrics/pbmc_modularity.csv', index_col=[0])
pbmc_sil = pd.read_csv('Data/cluster_metrics/pbmc_silhouette_scores.csv', index_col=[0])

sankey_dict = CellLayers(pbmc_exp, pbmc_meta, 
                         genes=['MS4A1', 'GNLY', 'CD3E', 'FCER1A'],
                         tri_coexpressed_genes=[['MS4A1', 'FCER1A', 'FCGR3A']])
```

### Seurat Implementation
```R
library(Seurat)
```

### Scanpy Implementation
```Python
import scanpy as sc 
```

## Documentation
Please consider citing Cell Layers if you used the application for your analysis.
