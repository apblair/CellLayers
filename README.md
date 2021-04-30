# Cell Layers: Uncovering clustering structure and knowledge in unsupervised single-cell transcriptomic analysis
#### Andrew P. Blair, Robert K. Hu, Katherine S. Pollard, Pawel F. Przytycki*, Irfan S. Kathiriya*, Benoit G. Bruneau*

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
$ docker pull bruneaulab/cell-layers:0.1
```

## Running the Container

To run the Cell Layers JupyterLab container please install Docker, open a terminal, and run the following command. 
```bash
$ docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work docker push bruneaulab/cell-layers:0.1
```

The Docker container was also converted to a Singularity container. 
```bash
$ singularity build cell-layers.sif docker://bruneaulab/cell-layers:0.1
$ singularity exec cell-layers.sif start.sh jupyter lab --port=9595
```


## Usage
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


## Documentation
Please consider citing Cell Layers if you used the application for your analysis.

## Authors Contributions
A.P.B. conceived and initiated the project. R.H assisted in analysis. B.G.B., P.F.P, and I.S.K. supervised A.P.B. I.S.K. provided datasets. K.S.P. advised. All authors commented on the manuscript.

## Acknowledgments
We thank the Cytoscape and scNetViz developers Alex Pico and Scooter Morris for their input on Plotly, Dan Carlin for his input on multi-resolution analysis, members of the CIRM Heart of Cells group, Gladstone Bioinformatics core, and Bruneau lab for discussions and comments. 

## Funding
California Institute for Regenerative Medicine (RB4-05901 to B.G.B)


