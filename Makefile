#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	make
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}


all: setup download_scripts orthofuser blast spades trinity shannon seqtk busco

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
ifeq "$(shell basename $(shell which orthofuser.py))" "orthofuser.py"
	@echo "ORTHOFUSER is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/macmanes-lab/OrthoFinder.git && export PATH=${DIR}/software/OrthoFinder/orthofinder:$$PATH
endif

blast:
ifeq "$(shell basename $(shell which blastp))" "blastp"
	@echo "BLASTP is already installed"
else
	@echo "blastp is not installed, installing now..."
	cd ${DIR}/software &$ curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz && tar -zxf ncbi-blast-2.6.0+-x64-linux.tar.gz
	export PATH=${DIR}/software/ncbi-blast-2.6.0+/bin:$$PATH
endif

spades:
ifeq "$(shell basename $(shell which spades.py))" "spades.py"
	@echo "SPAdes is already installed"
else
	cd ${DIR}/software && \
	curl -LO http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz && tar -zxf SPAdes-3.10.1-Linux.tar.gz && \
	export PATH=${DIR}/software/SPAdes-3.10.1-Linux/bin:$$PATH
endif

trinity:
ifeq "$(shell basename $(shell which Trinity))" "Trinity"
	@echo "Trinity is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/trinityrnaseq/trinityrnaseq.git && cd trinityrnaseq && make
	export PATH=$$PATH:${DIR}/software/trinityrnaseq
endif

shannon:
ifeq "$(shell basename $(shell which shannon.py))" "shannon.py"
	@echo "Shannon is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/sreeramkannan/Shannon.git
	export PATH=$$PATH:${DIR}/software/Shannon
endif

seqtk:
ifeq "$(shell basename $(shell which seqtk))" "seqtk"
	@echo "seqtk is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/lh3/seqtk.git && cd seqtk && make
	export PATH=$$PATH:${DIR}/software/seqtk
endif

busco:
ifeq "$(shell basename $(shell which run_BUSCO.py))" "run_BUSCO.py"
	@echo "BUSCO is already installed"
else
	cd ${DIR}/software && \
	git clone https://gitlab.com/ezlab/busco.git && cd busco && python setup.py install --user --prefix=
	export PATH=$$PATH:${DIR}/software/busco
endif
