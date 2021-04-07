FROM ubuntu:20.04
MAINTAINER Matthew MacManes (Matthew.MacManes@unh.edu)

#########
### Working dir
#########
RUN	mkdir build
WORKDIR /build

#########
### Setup the environment
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
    && locale-gen en_US.UTF-8 \
    && apt update && apt install -y build-essential git cmake python3 autoconf libbz2-dev liblzma-dev libxml2-dev libz-dev curl wget sudo nano \
    && apt -y autoremove


#########
### Create user orp
#########
RUN useradd -r -s /bin/bash -U -m -d /home/orp -p '' orp

#########
### Setup the user environment
########
ENV HOME /home/orp
RUN usermod -aG sudo orp \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#########
### Oyster River Protocol (ORP)
#########
WORKDIR /home/orp
USER orp
RUN git clone https://github.com/macmanes-lab/Oyster_River_Protocol.git && cd Oyster_River_Protocol && sudo make && sed -i  "s_ubuntu_orp_g" /home/orp/Oyster_River_Protocol/software/config.ini

#########
### Path
#########
RUN echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/anaconda/install/bin' >> /home/orp/.profile \
    && echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/OrthoFinder/orthofinder' >> /home/orp/.profile \
    && echo 'export TRINITY_HOME=/home/orp/Oyster_River_Protocol/software/trinityrnaseq-v2.12.0/' >> /home/orp/.profile \
    && echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/SPAdes-3.15.2-Linux/bin/' >> /home/orp/.profile \
    && echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/orp-transrate' >> /home/orp/.profile \
    && echo 'PATH=$PATH:/home/orp/Oyster_River_Protocol/software/transabyss' >> /home/orp/.profile \
    && /bin/bash -c "source /home/orp/.profile"
ENV PATH="/home/orp/Oyster_River_Protocol/software/transabyss:/home/orp/Oyster_River_Protocol/software/anaconda/install/bin:/home/orp/Oyster_River_Protocol/software/trinityrnaseq-v2.12.0/:/home/orp/Oyster_River_Protocol/software/orp-transrate:/home/orp/Oyster_River_Protocol/software/SPAdes-3.15.2-Linux/bin/:/home/orp/Oyster_River_Protocol/software/OrthoFinder/orthofinder:/home/orp/Oyster_River_Protocol/software/anaconda/install/bin:${PATH}"

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
USER root
RUN chown -R orp:orp  /home/orp/Oyster_River_Protocol/busco_dbs/
USER orp
RUN /bin/bash -c "source /home/orp/.profile" \
    && /bin/bash -c "source /home/orp/Oyster_River_Protocol/software/anaconda/install/etc/profile.d/conda.sh"
