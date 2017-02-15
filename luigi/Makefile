.PHONY: run containers

all: run

containers:
	$(MAKE) -C containers


first_pass_meerkat_ar1.rfis:
	tar zxvf first_pass_meerkat_ar1.rfis.tgz

rfi_mask.pickle:
	tar zxvf rfi_mask.pickle.tgz

run: containers first_pass_meerkat_ar1.rfis rfi_mask.pickle
	PYTHONPATH='.' luigi --module vermeertasks --logging-conf-file logging.conf WscleanTask  --local-scheduler
