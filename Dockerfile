FROM ubuntu
# docker build -t 10.168.2.67:5000/perbase-yy:1.0 .
# docker push 10.168.2.67:5000/perbase-yy:1.0
# Maintainer: wangxz
# Date: 2021-03-02

### create required directory
RUN mkdir -p /yunying/software/

RUN groupadd -g 10001 yunying \
    && useradd -g yunying -u 10001 yunying \
    && chmod 775 /yunying

COPY target/perbase /usr/local/bin/

COPY script/*.py /opt/

USER yunying

ENTRYPOINT ["/bin/bash", "-c"]
