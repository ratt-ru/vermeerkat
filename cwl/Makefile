.PHONY: run images

all: run


.virtualenv:
	virtualenv -p python2 .virtualenv

requirements: .virtualenv
	.virtualenv/bin/pip install -r requirements.txt

image:
	$(MAKE) -C images

run: images requirements
	.virtualenv/bin/cwltoil workflow.cwl job.cwl
