VERSION=0.1
DOCKER_REPO=vermeerkat/rfimasker:$(VERSION)


.PHONY: build clean run


all: kliko

unpack:
	tar zxvf rfi_mask.pickle.tgz

build: unpack
	docker build -t ${DOCKER_REPO} .

force-build:
	docker build --no-cache -t ${DOCKER_REPO} .

clean:
	docker rmi ${DOCKER_REPO}


run: build
	docker run -ti ${DOCKER_REPO}

kliko: build 
	kliko-run ${DOCKER_REPO} --mask rfi_mask.pickle
