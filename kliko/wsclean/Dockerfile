FROM kernsuite/base:dev
RUN docker-apt-install wsclean
RUN pip install kliko==0.7.1
ADD kliko.yml /
ADD kliko /
CMD /usr/bin/wsclean
