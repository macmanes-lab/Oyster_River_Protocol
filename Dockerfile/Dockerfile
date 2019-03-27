FROM delhomme/upscb-ngs-preprocessing
MAINTAINER Nicolas Delhomme (nicolas.delhomme@umu.se)

#########
### Working dir
#########
RUN	mkdir build
WORKDIR /build

#########
### Aptitude packages
#########
RUN apt update && apt install -y automake autotools-dev cmake emboss hmmer infernal libboost-all-dev \
libfontconfig1 libsm6 liburi-escape-xs-perl liburi-perl libxrender1 mcl ncbi-blast+ parallel python-numpy \
ruby unzip && gem install crb-blast


#########
### Python modules
#########
RUN pip install --upgrade pip


#########
### Oyster River Protocol (ORP)
#########
WORKDIR /home/training
RUN git clone https://github.com/macmanes-lab/Oyster_River_Protocol.git
RUN cd Oyster_River_Protocol && make
USER root

#########
### Path
#########
RUN echo 'PATH=$PATH:/home/training/miniconda/bin' >> /home/training/.profile
RUN echo 'PATH=$PATH:/home/training/miniconda/bin' >> /home/training/.bash_profile
ENV PATH=/home/training/miniconda/bin:${PATH}


#########
### Clean up
#########
WORKDIR /
RUN rm -rf /build

#########
### Ports and CMD
#########
CMD ["/usr/bin/supervisord","-c","/etc/supervisor/conf.d/supervisord.conf"]
