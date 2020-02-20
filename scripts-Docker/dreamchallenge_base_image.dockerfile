FROM ubuntu:latest
MAINTAINER Salem Malikic <salemmalikic05@gmail.com>
USER root
RUN apt-get -y update
RUN apt-get install -y python python-dev python-distribute python-pip
RUN apt-get -y install r-base r-base-dev
RUN Rscript -e 'install.packages("ggplot2", repo="http://cran.ma.imperial.ac.uk/")'
RUN Rscript -e 'install.packages("RColorBrewer", repo="http://cran.ma.imperial.ac.uk/")'
RUN Rscript -e 'install.packages("DPpackage", repo="http://cran.ma.imperial.ac.uk/")'
RUN Rscript -e 'install.packages("lpSolve", repo="http://cran.ma.imperial.ac.uk/")'
RUN apt-get -y install git
RUN mkdir /opt/CTPsingle
