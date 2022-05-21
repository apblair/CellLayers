FROM jupyter/datascience-notebook:latest

USER root
RUN conda install -c conda-forge r-seurat
RUN R -e "install.packages(c('cluster','enrichR'), repos='http://cran.us.r-project.org')"
# install.packages('BiocManager')
# BiocManager::install('limma')
USER $NB_UID
# RUN pip install CellLayers
