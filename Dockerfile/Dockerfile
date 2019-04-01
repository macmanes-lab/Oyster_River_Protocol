FROM ubuntu:18.04
MAINTAINER Matthew MacManes (Matthew.MacManes@unh.edu)

#########
### Working dir
#########
RUN	mkdir build
WORKDIR /build

#########
### Setup the envirenment
#########
ENV DEBIAN_FRONTEND noninteractive
ENV LC_ALL en_GB.UTF-8
ENV LANG en_GB.UTF-8
ENV LANGUAGE en_GB:en

#########
### Aptitude packages
#########
RUN apt update \
    && apt -y upgrade \
    && apt install -y --reinstall language-pack-en \
    && locale-gen en_GB.UTF-8
RUN apt update && apt install -y build-essential git python python-pip libxml2-dev libz-dev curl wget


#########
### Python modules
#########
RUN pip install --upgrade pip

#########
### Create user training
#########
RUN useradd -r -s /bin/bash -U -m -d /home/training -p '' training

#########
### Setup the user environment
########
ENV HOME /home/training
RUN usermod -aG sudo,audio,video training \
&& echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#########
### Oyster River Protocol (ORP)
#########
WORKDIR /home/training
USER training
RUN git clone https://github.com/macmanes-lab/Oyster_River_Protocol.git
RUN cd Oyster_River_Protocol && make
RUN sed -i  "s_ubuntu_training_g" /home/training/Oyster_River_Protocol/software/config.ini

#########
### Path
#########
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/anaconda/install/bin' >> /home/training/.profile
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/OrthoFinder/orthofinder' >> /home/training/.profile
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/orp-transrate' >> /home/training/.profile
RUN echo 'PATH=$PATH:/home/training/Oyster_River_Protocol/software/transabyss' >> /home/training/.profile
RUN /bin/bash -c "source /home/training/.profile"
ENV PATH="/home/training/Oyster_River_Protocol/software/transabyss:/home/training/Oyster_River_Protocol/software/orp-transrate:/home/training/Oyster_River_Protocol/software/OrthoFinder/orthofinder:/home/training/Oyster_River_Protocol/software/anaconda/install/bin:${PATH}"

#########
### Clean up
#########
USER root
WORKDIR /
RUN rm -rf /build

#########
### Setup init workspace
#########
WORKDIR $HOME
USER training
