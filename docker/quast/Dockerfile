FROM tpesout/hpp_base:latest
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu

ENV PATH="/opt/python/Python-3.6.0/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y pkg-config libfreetype6-dev libpng-dev python-matplotlib python-simplejson libjpeg-dev zlib1g-dev


### QUAST
WORKDIR /opt/quast
RUN wget https://github.com/ablab/quast/releases/download/quast_5.3.0/quast-5.3.0.tar.gz && \
    tar -xzf quast-5.3.0.tar.gz && \
    rm quast-5.3.0.tar.gz && \
    cd quast-5.3.0 && \
    pip3 install --upgrade pip && \
    pip3 install --upgrade setuptools==59.6.0 matplotlib==3.3.4 simplejson==3.17.6 joblib==0.14.1 pyparsing==2.4.7 pysam==0.19.1 && \
    ./setup.py install && \
    ./quast.py -h

# locale configuration (there is an issue in python3 parsing a reference sequence)
RUN apt-get install --reinstall -y locales && \
    locale-gen en_US.UTF-8 

WORKDIR /data

