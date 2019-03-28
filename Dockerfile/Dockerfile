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
RUN apt update && apt install -y build-essential git python python-pip libxml2-dev libz-dev


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
RUN sed -i  "s_ubuntu_$(whoami)_g" $HOME/Oyster_River_Protocol/software/config.ini
USER root

#########
### Path
#########
RUN cat pathfile >> ~/.profile
RUN cat pathfile >> ~/.bash_profile
RUN source ~/.profile


#########
### Clean up
#########
WORKDIR /
RUN rm -rf /build

#########
### Ports and CMD
#########
CMD ["/usr/bin/supervisord","-c","/etc/supervisor/conf.d/supervisord.conf"]
