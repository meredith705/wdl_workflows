FROM ubuntu:20.04
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu

RUN mkdir -p /home/apps

RUN cd /home/apps && \
    apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget make build-essential cmake \
    protobuf-compiler pkg-config libprotobuf-dev libjansson-dev libhts-dev


RUN cd /home/apps && \
    git clone https://github.com/rlorigro/GFAse.git && \
    cd GFAse && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 8 
ENV PATH="/home/apps/GFAse/build:${PATH}"

RUN mkdir -p /data
WORKDIR /data
