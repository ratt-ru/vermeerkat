VERSION=0.1
DOCKER_REPO=vermeerkat/casasouthcal:$(VERSION)


.PHONY: build clean run build-patch


all: build


build-patch:
	docker build -t casasouthcaltmp -f patch.docker .

fluxcalibrator.data/: build-patch
	docker run casasouthcaltmp tar -c -C /casa-release-4.7.0-el6/data/nrao/VLA/standards fluxcalibrator.data | tar x

PerleyButler2013Coeffs/: build-patch
	docker run casasouthcaltmp tar -c -C /casa-release-4.7.0-el6/data/nrao/VLA/standards PerleyButler2013Coeffs | tar x

build: fluxcalibrator.data PerleyButler2013Coeffs
	docker build -t ${DOCKER_REPO} .

force-build:
	docker build --no-cache -t ${DOCKER_REPO} .

clean:
	docker rmi ${DOCKER_REPO}


run: build
	docker run -ti ${DOCKER_REPO} bash
