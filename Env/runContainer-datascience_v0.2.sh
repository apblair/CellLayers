docker run --rm -p 10000:8888 \
	--user root \
	-e JUPYTER_ENABLE_LAB=yes \
	-e CHOWN_HOME=yes \
	-v "$PWD":/home/jovyan/work apblair/cell-layers:v0.2