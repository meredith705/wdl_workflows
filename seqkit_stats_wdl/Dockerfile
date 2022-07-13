FROM ubuntu:20.04
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu


RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget build-essential python3 python3-dev python3-pip python3-matplotlib tabix

RUN pip3 install numpy matplotlib cython==0.29

WORKDIR /opt/

RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.0.0/seqkit_linux_amd64.tar.gz && tar zxvf seqkit_linux_amd64.tar.gz \
&& mv seqkit /usr/local/bin/ && rm seqkit_linux_amd64.tar.gz


WORKDIR /data