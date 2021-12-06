docker run -it -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes apblair/cell-layers:v0.2
# -e JUPYTER_ENABLE_LAB=yes 
# Inside the container run 
# $ mamba install -c conda-forge jupyterlab
# $ jupyter lab --ip 0.0.0.0 --no-browser --allow-root
