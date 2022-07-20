FROM jupyter/datascience-notebook:latest

COPY install-SetupCellLayers.R /tmp/install-SetupCellLayers.R

USER root
RUN conda install -c conda-forge r-seurat
RUN R -e "install.packages(c('cluster','enrichR','BiocManager'), repos='http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('limma')"
RUN Rscript /tmp/install-SetupCellLayers.R

USER $NB_UID
RUN pip install --user CellLayers==0.2.2