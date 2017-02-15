FROM kernsuite/base:1
RUN pip install kliko==0.8
RUN docker-apt-install wsclean
ADD kliko.yml /
ADD kliko /
CMD /usr/bin/wsclean
