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
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en

#########
### Aptitude packages
#########
RUN apt update \
    && apt -y upgrade \
    && apt install -y --reinstall language-pack-en \
    && locale-gen en_US.UTF-8
RUN apt update && apt install -y build-essential git python python-pip libxml2-dev libz-dev curl wget sudo nano


#########
### Python modules
#########
RUN pip install --upgrade pip

#########
### Create user orp
#########
RUN useradd -r -s /bin/bash -U -m -d /home/orp -p '' orp

#########
### Setup the user environment
########
ENV HOME /home/orp
RUN usermod -aG sudo,audio,video orp \
&& echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#########
### Oyster River Protocol (ORP)
#########
WORKDIR /home/orp
USER orp
RUN git clone https://github.com/macmanes-lab/Oyster_River_Protocol.git
RUN cd Oyster_River_Protocol && sudo make
RUN sed -i  "s_ubuntu_orp_g" /home/orp/Oyster_River_Protocol/software/config.ini

#########
### Path
#########
RUN echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/anaconda/install/bin' >> /home/orp/.profile
RUN echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/OrthoFinder/orthofinder' >> /home/orp/.profile
RUN echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/orp-transrate' >> /home/orp/.profile
RUN echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/transabyss' >> /home/orp/.profile
RUN /bin/bash -c "source /home/orp/.profile"
ENV PATH="/home/orp/Oyster_River_Protocol/software/transabyss:/home/orp/Oyster_River_Protocol/software/orp-transrate:/home/orp/Oyster_River_Protocol/software/OrthoFinder/orthofinder:/home/orp/Oyster_River_Protocol/software/anaconda/install/bin:${PATH}"

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
USER orp
RUN /bin/bash -c "source activate orp"
