VERSION=0.1
DOCKER_REPO=vermeerkat/downobs:$(VERSION)


.PHONY: build clean run kliko


all: kliko

force-build:
	docker build --no-cache -t ${DOCKER_REPO} .

build:
	docker build -t ${DOCKER_REPO} .


clean:
	docker rmi ${DOCKER_REPO}


run: build
	docker run -ti ${DOCKER_REPO}

kliko: build
	kliko-run ${DOCKER_REPO} --url http://kat-archive.kat.ac.za/archive/data/comm/2012/09/29/1348904076.h5 --filename 1348904076.h5
