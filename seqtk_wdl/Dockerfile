FROM ubuntu:20.04


MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu


RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev zlib1g-dev && \
	apt-get clean





# seqtk clone and install
RUN git clone https://github.com/lh3/seqtk.git /opt/seqtk && \
	cd /opt/seqtk && \
	make && \
	cp seqtk /usr/bin





