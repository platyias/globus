
FROM ubuntu:trusty
MAINTAINER German Plata <platyias@yahoo.com>
WORKDIR  /data
ADD . /data
ENV CLASSPATH /data/scripts
RUN apt-get update && \
apt-get -y --no-install-recommends install g++ default-jdk make nano libpython2.7 python2.7 python-pip python-dev && \
cd /data/CGibbs && \
make && \
cd /data/scripts && \
javac GeneClusteringSplitTL.java && \
pip install numpy pandas && \
cd /data && \
mkdir diamond && \
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz && \
tar -xvf diamond-linux64.tar.gz -C diamond && \



