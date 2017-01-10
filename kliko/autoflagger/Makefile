VERSION=0.1
DOCKER_REPO=vermeerkat/autoflagger:$(VERSION)


.PHONY: build clean run


all: kliko



build:
	docker build -t ${DOCKER_REPO} .

force-build:
	docker build --no-cache -t ${DOCKER_REPO} .

clean:
	docker rmi ${DOCKER_REPO}


run: build
	docker run -ti ${DOCKER_REPO}

first_pass_meerkat_ar1.rfis:
	tar zxvf first_pass_meerkat_ar1.rfis.tgz

kliko: build first_pass_meerkat_ar1.rfis
	kliko-run ${DOCKER_REPO} --strategy first_pass_meerkat_ar1.rfis
