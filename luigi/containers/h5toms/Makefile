VERSION=0.1
DOCKER_REPO=vermeerkat/h5toms:$(VERSION)


.PHONY: build clean run


all: kliko

input/1348904076.h5:
	curl http://kat-archive.kat.ac.za/archive/data/comm/2012/09/29/1348904076.h5 --output input/1348904076.h5

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

kliko: build input/1348904076.h5
	kliko-run ${DOCKER_REPO} --pattern 1348904076.h5
