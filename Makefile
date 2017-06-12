#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	make
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}


prep: setup run_scripts

.DELETE_ON_ERROR:
.PHONY:report

setup:
	mkdir -p ${DIR}/scripts
	mkdir -p ${DIR}/software
	mkdir -p ${DIR}/shared


download_scripts:
	cd ${DIR}/scripts && \
	curl -LO https://raw.githubusercontent.com/macmanes-lab/general/master/filter.py

orthofuser:
	cd ${DIR}/software && \
	git clone https://github.com/macmanes-lab/OrthoFinder.git && export PATH=${DIR}/software/OrthoFinder/orthofinder:$PATH

blast:
	ifndef blastp && \
		cd ${DIR}/software && \
		curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz && tar -zxf ncbi-blast-2.6.0+-x64-linux.tar.gz && \
		export PATH=${DIR}/software/ncbi-blast-2.6.0+/bin:$PATH
spades:
	ifndef spades.py && \
	curl -LO http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz && tar -zxf SPAdes-3.10.1-Linux.tar.gz && \
	export PATH=${DIR}/software/SPAdes-3.10.1-Linux/bin:$PATH
