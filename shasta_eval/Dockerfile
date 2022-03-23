FROM ubuntu:20.04
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu

RUN mkdir -p /home/apps

RUN cd /home/apps && \
    apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget make build-essential cmake \
    protobuf-compiler pkg-config libprotobuf-dev libjansson-dev libhts-dev

ADD https://api.github.com/repos/rlorigro/GFAse/git/refs/heads/wdl_debug version.json

RUN cd /home/apps && \
    git clone -b wdl_debug https://github.com/rlorigro/GFAse.git && \
    cd GFAse && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 8 
    
ENV PATH="/home/apps/GFAse/build:${PATH}"
ENV PATH="/home/apps/GFAse/data:${PATH}"

RUN mkdir -p /data
WORKDIR /data
