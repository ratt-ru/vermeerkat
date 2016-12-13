VERSION=0.1
DOCKER_REPO=kernsuite/h5toms:$(VERSION)


.PHONY: build clean run


all: kliko


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

kliko: build
	kliko-run ${DOCKER_REPO}
