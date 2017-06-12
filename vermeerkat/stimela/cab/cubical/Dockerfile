FROM stimela/base:0.2.8
MAINTAINER <sphemakh@gmail.com>
ENV TERM xterm
ADD src /code
ENV LOGFILE ${OUTPUT}/logfile.txt
ENTRYPOINT ["/bin/bash"]
CMD ["-c", "python /code/run.py 2>&1 | tee -a $LOGFILE ; (exit ${PIPESTATUS})"]