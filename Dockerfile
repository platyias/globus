
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
pip install numpy && \
pip install pandas
