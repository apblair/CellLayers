API
=====

.. _SetupCellLayers:

SetupCellLayers
---------------
To provide options for generating the input data to Cell Layers, an R library (SetupCellLayers) is available to generate a cell-by-resolution parameter matrix from a scRNA-seq kNN graph using the popular Seurat SLM and SNN clustering methods.

===================================  ====================  
Function                             Summary            
===================================  ====================
compute_multi_res_search             Compute clustering for all user defined resolution parameters.   
compute_silhouette_scores            Compute silhouette scores using the PCA embedding space for each resolution's communities. 
compute_enrichment                   Compute an enrichment analysis of each community using their differentially expressed genes for each resolution parameter.
===================================  ====================

.. _CellLayers:

CellLayers
----------
===================================  ====================
Function                             Summary
===================================  ====================
build_sankey
build_enrichment_sankey
build_coexpression_sankey
===================================  ====================
