DOCKER_REPO=vermeerkat/h5toms:latest


.PHONY: build clean run


all: run

unzip:
	tar zxvf blank.ms.tgz

build: unzip
	docker build -t ${DOCKER_REPO} .

force-build:
	docker build --no-cache -t ${DOCKER_REPO} .

clean:
	docker rmi ${DOCKER_REPO}

run: build
	docker run -ti ${DOCKER_REPO}