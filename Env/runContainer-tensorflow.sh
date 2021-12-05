docker run --rm -p 10000:8888 \
	-e JUPYTER_ENABLE_LAB=yes \
	-e CHOWN_HOME=yes \
	-e CHOWN_HOME_OPTS='R' \
	-v "$PWD":/home/jovyan/work docker.synapse.org/syn26380219/cell-layers:v0.1
