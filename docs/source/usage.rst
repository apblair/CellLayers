Usage
=====

.. _installation:

Installation
------------

The **CellLayers** Python module can be installed via pip. Please note, the Python module requires the cluster metric and enrichment data to be generated independently or via our R library **SetupCellLayers**. 

.. code-block:: console

   (.venv) $ pip install CellLayers

The **SetupCellLayers** R library can be installed via devtools.

.. code-block:: R

   > library(devtools)
   > options(timeout=9999999)
   > devtools::install_github("apblair/CellLayers/SetupCellLayers")
   > library(SetupCellLayers)