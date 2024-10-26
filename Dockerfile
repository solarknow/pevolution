FROM ubuntu:latest
MAINTAINER Omniscient Sun

RUN echo "0.5" > /version

ADD . /

RUN ./pevolution-install.sh
