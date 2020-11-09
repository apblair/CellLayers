# Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis
Andrew P. Blair, Robert K. Hu, Irfan S. Kathiriya, Katherine S. Pollard, Pawel F. Przytycki, Benoit G. Bruneau

## Motivation
Unsupervised clustering of single-cell transcriptomics is a powerful method for identifying cell populations. Static visualization techniques for single-cell clustering only display results for a single resolution parameter. Analysts will often evaluate more than one resolution parameter, but then only report one.

## Results
Cell Layers is an interactive Sankey tool for the quantitative investigation of gene expression, coexpression, biological processes, and cluster integrity across clustering resolutions. Cell Layers enhances the interpretability of single-cell clustering by linking molecular data and cluster evaluation metrics to provide novel insight into cell populations.


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

exp_df = pd.read_csv()
meta_df = pd.read_csv()
```

### Seurat implementation
```R
library(Seurat)
```

### Scanpy implementation
```Python
import 
```

## Documentation
Please consider citing Cell Layers if you used the application for your analysis.
