Usage
=====

.. _installation:

Installation
------------

The **CellLayers** Python module can be installed via pip. Please note, the Python module requires the cluster metric and enrichment data to be generated independently or via our R library **SetupCellLayers**. 

.. code-block:: console

   $ pip install CellLayers

The **SetupCellLayers** R library can be installed via devtools.

.. code-block:: R

   > library(devtools)
   > options(timeout=9999999)
   > devtools::install_github("apblair/CellLayers/SetupCellLayers")
   > library(SetupCellLayers)

Both CellLayers and SetupCellLayers are fully containerized via Docker and Singularity, which are extended from the jupyter/datascience-notebook image.

.. code-block:: console
   
   $ docker pull apblair/cell-layers:v0.2
   $ docker run -it --rm -p 10000:8888 -v "${PWD}":/home/jovyan/work apblair/cell-layers:v0.2

.. code-block:: console
   
   $ singularity pull --arch amd64 library://apblair/single-cell-tools/cell-layers:v0-2
   $ singularity exec cell-layers_v0-2.sif start.sh jupyter lab --port=9595

.. toctree::
   :maxdepth: 2

   notebooks/PBMC_SetupCellLayer.ipynb