FROM ubuntu:18.04
MAINTAINER Matthew MacManes (Matthew.MacManes@unh.edu)

#########
### Working dir
#########
RUN	mkdir build
WORKDIR /build

#########
### Aptitude packages
#########
RUN apt update && apt install -y build-essential git python python-pip libxml2-dev libz-dev curl wget


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
RUN sed -i  "s_ubuntu_training_g" /home/training/Oyster_River_Protocol/software/config.ini
USER root

#########
### Path
#########
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/anaconda/install/bin' >> /home/training/.profile
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/OrthoFinder/orthofinder' >> /home/training/.profile
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/orp-transrate' >> /home/training/.profile
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/transabyss' >> /home/training/.profile
ENV PATH="/home/training/Oyster_River_Protocol/software/transabyss:/home/training/Oyster_River_Protocol/software/orp-transrate:/home/training/Oyster_River_Protocol/software/OrthoFinder/orthofinder:/home/training/Oyster_River_Protocol/software/anaconda/install/bin:${PATH}"
RUN /bin/bash -c "source /home/training/.profile"
USER root

#########
### Clean up
#########
WORKDIR /
RUN rm -rf /build
