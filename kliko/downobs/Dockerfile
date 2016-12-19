FROM kernsuite/base:dev
RUN docker-apt-install curl
RUN pip install kliko==0.7.1
ADD kliko /
ADD kliko.yml /
CMD /usr/bin/curl