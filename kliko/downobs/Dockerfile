FROM kernsuite/base:1
RUN pip install kliko==0.8
RUN docker-apt-install curl
ADD kliko /
ADD kliko.yml /
CMD /usr/bin/curl
