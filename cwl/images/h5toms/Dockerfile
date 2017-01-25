FROM kernsuite/base:1
RUN docker-apt-install python-katdal python-casacore
ADD blank.ms /var/kat/static/blank.ms
CMD /usr/bin/h5toms.py
