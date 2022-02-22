FROM ubuntu:20.04
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu


RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget make build-essential cmake \
    protobuf-compiler pkg-config libprotobuf-dev libjansson-dev libhts-dev


WORKDIR /opt/
RUN git clone https://github.com/rlorigro/GFAse.git && \
    cd GFAse && \
    cmake . && \
    make -j 8 


WORKDIR /data
