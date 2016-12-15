DOCKER_REPO=vermeerkat/wsclean:0.1


.PHONY: build force-build clean run


all: kliko


build:
	docker build -t ${DOCKER_REPO} .

force-build:
	docker build --no-cache -t ${DOCKER_REPO} .

clean:
	docker rmi ${DOCKER_REPO}

run: build
	docker run -ti ${DOCKER_REPO} bash

kliko: build
	kliko-run ${DOCKER_REPO} 
