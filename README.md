# Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis

Cell Layers is an interactive Sankey tool for the quantitative investigation of gene expression, coexpression, biological processes, and cluster integrity across clustering resolutions.

## Installation
Cell Layers may be installed via pip after pulling the repo.

```bash
$ git clone git@github.com:apblair/CellLayers.git
$ cd CellLayers
$ pip install .
```

```bash
$ docker pull bruneaulab/cell-layers:0.1
```

## Running the Container

```bash
$ docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work bruneaulab/cell-layers:0.1
```

For HPC clusters that support Linux containers, the Docker container may be pulled and built as a Singularity container. After building the Singularity container, an interactive session can be run with these commands:
```bash
$ singularity build cell-layers.sif docker://bruneaulab/cell-layers:0.1
$ singularity exec cell-layers.sif start.sh jupyter lab --port=9595
```

## Usage

Please unzip the tutorial's expression data.

```bash
$ cd CellLayers/Data
$ unzip PBMC_exp.csv.zip 
```

Next, open CellLayers/Notebooks/PBCM_Tutorial.ipynb in a Jupyter environment and run the following cell:

```python
import CellLayers
import pandas as pd

pbmc_exp = pd.read_csv('Data/PBMC_exp.csv', index_col=[0]) # cell by gene expression matrix
pbmc_meta = pd.read_csv('Data/PBMC_meta.csv', index_col=[0]) # cell by resolution matrix
pbmc_modularity = pd.read_csv('Data/pbmc_modularity.csv', index_col=[0])
pbmc_sil = pd.read_csv('Data/pbmc_silhouette_scores.csv', index_col=[0])

sankey = CellLayers.run(exp_df, meta_df, modularity=mod_df, silhouette=sil_df, genes=['CD3E'])
sankey.show()
```
<image src="Images/example.png">
    
**Fig. 1 (A)** PBMC mulit-resolution analysis from 0.1 to 0.5. Edges are painted by CD3E, which is a marker gene for CD8+ T, Memory CD4 T, and Naive CD4 T cells. Nodes are painted by Silhouette score. The lower Sihouette values indicate samples are near the decision boundary of neighboring clusters.
    
**Fig. 1 B**

**Fig. 1 C**
  
    
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


