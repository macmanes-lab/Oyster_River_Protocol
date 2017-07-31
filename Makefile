#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	make
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
rcorrPath = `which rcorrector`

all: setup brew orthofuser rcorrector blast spades trinity shannon seqtk busco trimmomatic transrate postscript

.DELETE_ON_ERROR:
.PHONY:report

setup:
	@mkdir -p ${DIR}/scripts
	@mkdir -p ${DIR}/software
	@mkdir -p ${DIR}/shared
	@rm -f pathfile

brew:
ifeq "$(shell basename $(shell which brew))" "brew"
	@echo "BREW is already installed"
else
	$error("*** BREW MUST BE PROPERLY INSTALLED BEFORE YOU CAN PROCEED, SEE: http://angus.readthedocs.io/en/2016/linuxbrew_install.html ***")
endif

rcorrector:
	@if [ $$(basename $$(which rcorrector)) == 'rcorrector' ];\
	then\
		@echo "Rcorrector is already installed"
	else\
		brew install rcorrector
	fi

orthofuser:
ifeq "$(shell basename $(shell which orthofuser.py))" "orthofuser.py"
	@touch pathfile
	@echo "ORTHOFUSER is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/macmanes-lab/OrthoFinder.git
	@echo ${DIR}/software/OrthoFinder/orthofinder | tee -a pathfile
endif

blast:
ifeq "$(shell basename $(shell which blastp))" "blastp"
	@echo "BLASTP is already installed"
else
	@echo "blastp is not installed, installing now..."
	cd ${DIR}/software &$ curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz && tar -zxf ncbi-blast-2.6.0+-x64-linux.tar.gz
	@echo ${DIR}/software/ncbi-blast-2.6.0+/bin | tee -a pathfile
endif

spades:
ifeq "$(shell basename $(shell which spades.py))" "spades.py"
	@echo "SPAdes is already installed"
else
	cd ${DIR}/software && \
	curl -LO http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz && tar -zxf SPAdes-3.10.1-Linux.tar.gz && \
	@echo ${DIR}/software/SPAdes-3.10.1-Linux/bin | tee -a pathfile
endif

trinity:
ifeq "$(shell basename $(shell which Trinity))" "Trinity"
	@echo "TRINITY is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/trinityrnaseq/trinityrnaseq.git && cd trinityrnaseq && make
	@echo PATH=$$PATH:${DIR}/software/trinityrnaseq | tee -a pathfile
endif

shannon:
ifeq "$(shell basename $(shell which shannon.py))" "shannon.py"
	@echo "SHANNON is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/sreeramkannan/Shannon.git
	@echo PATH=$$PATH:${DIR}/software/Shannon | tee -a pathfile
endif

seqtk:
ifeq "$(shell basename $(shell which seqtk))" "seqtk"
	@echo "SEQTK is already installed"
else
	cd ${DIR}/software && \
	git clone https://github.com/lh3/seqtk.git && cd seqtk && make
	@echo PATH=$$PATH:${DIR}/software/seqtk | tee -a pathfile
endif

busco:
ifeq "$(shell basename $(shell which run_BUSCO.py))" "run_BUSCO.py"
	@echo "BUSCO is already installed"
else
	cd ${DIR}/software && \
	git clone https://gitlab.com/ezlab/busco.git && cd busco && python setup.py install --user --prefix=
	@echo PATH=$$PATH:${DIR}/software/busco/scripts | tee -a pathfile
endif

trimmomatic:
ifeq "$(shell basename $(shell which trimmomatic))" "trimmomatic"
	@echo "TRIMMOMATIC is already installed"
else
	brew install trimmomatic
endif

transrate:
ifeq "$(shell basename $(shell which transrate))" "transrate"
	@echo "TRANSRATE is already installed"
else
	cd ${DIR}/software && \
	curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz && tar -zxf transrate-1.0.3-linux-x86_64.tar.gz
	@echo PATH=$$PATH:${DIR}/software/transrate-1.0.3-linux-x86_64 | tee -a pathfile
endif

postscript:
	@printf "\n\n*** The following location(s), if any print, need to be added to your PATH ***\n\n"
	@cat pathfile
	@printf "\n\n\n"
