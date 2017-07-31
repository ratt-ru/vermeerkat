FROM stimela/base:0.2.9
MAINTAINER <sphemakh@gmail.com>
ENV TERM xterm
ADD src /code
ADD cubical /src/cubical

RUN cd /src/cubical && \
    git checkout jsk-devel && \
    pip install -r requirements.txt && \
    pip install .

ENV LOGFILE ${OUTPUT}/logfile.txt
ENTRYPOINT ["/bin/bash"]
CMD ["-c", "python /code/run.py 2>&1 | tee -a $LOGFILE ; (exit ${PIPESTATUS})"]