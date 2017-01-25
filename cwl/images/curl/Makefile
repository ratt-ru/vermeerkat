DOCKER_REPO=vermeerkat/curl:latest


.PHONY: build clean run


all: run


force-build:
	docker build --no-cache -t ${DOCKER_REPO} .


build:
	docker build -t ${DOCKER_REPO} .


clean:
	docker rmi ${DOCKER_REPO}


run: build
	docker run -ti ${DOCKER_REPO} curl -O -L -C - http://kat-archive.kat.ac.za/archive/data/comm/2012/09/29/1348904076.h5
