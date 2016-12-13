FROM kernsuite/base:dev
RUN docker-apt-install python-katdal python-casacore
RUN pip install kliko==0.7.1
ADD kliko /
ADD kliko.yml /
ADD blank.ms /var/kat/static/blank.ms
