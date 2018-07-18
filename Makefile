#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	make
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CONDAROOT = ${DIR}/software/anaconda/install/
orthopath := $(shell which ${DIR}/software/OrthoFinder/orthofinder/orthofuser.py 2>/dev/null)
orthufuserversion = $(shell orthofuser.py --help | grep "OrthoFinder version" | awk '{print $$3}')
transrate := $(shell which ${DIR}/software/orp-transrate 2>/dev/null)
conda := $(shell which ${DIR}/software/anaconda/install/bin/activate 2>/dev/null)


all: setup conda orthofuser transrate shmlast_data busco_data postscript

.DELETE_ON_ERROR:

setup:
	@mkdir -p ${DIR}/scripts
	@mkdir -p ${DIR}/shared
	@rm -f pathfile

conda:
ifdef conda
	@echo "conda is already installed"
else
	mkdir -p ${DIR}/software/anaconda
	cd ${DIR}/software/anaconda && curl -LO https://repo.anaconda.com/archive/Anaconda3-5.1.0-Linux-x86_64.sh
	cd ${DIR}/software/anaconda && bash Anaconda3-5.1.0-Linux-x86_64.sh -b -p ${DIR}/software/anaconda/install
	( \
       source ${DIR}/software/anaconda/install/bin/activate; \
       conda update -y -n base conda; \
			 conda install -y --file <(curl https://raw.githubusercontent.com/macmanes-lab/Oyster_River_Protocol/conda/environment.yml); \
			 source deactivate; \
  )
	mkdir -p ${DIR}/software/shmlast && cd ${DIR}/software/shmlast && curl -LO ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && gzip -d uniprot_sprot.fasta.gz
endif

shmlast_data:
ifdef shmlast_data
	@echo "shmlast_data is already installed"
else
	mkdir -p ${DIR}/software/shmlast && cd ${DIR}/software/shmlast && curl -LO ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && gzip -d uniprot_sprot.fasta.gz
endif

busco_data:
	mkdir $HOME/Oyster_River_Protocol/busco_dbs && cd $HOME/Oyster_River_Protocol/busco_dbs
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz && tar -zxf eukaryota_odb9.tar.gz
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/arthropoda_odb9.tar.gz
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/insecta_odb9.tar.gz
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/tetrapoda_odb9.tar.gz
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/aves_odb9.tar.gz
	cd $HOME/Oyster_River_Protocol/busco_dbs && wget http://busco.ezlab.org/v2/datasets/mammalia_odb9.tar.gz

transrate:brew
ifdef transrate
	@echo "transrate is already installed"
else
	cd ${DIR}/software && tar -zxf orp-transrate.tar.gz
	@echo PATH=\$$PATH:${DIR}/software/orp-transrate >> pathfile
endif

orthofuser:
ifdef orthopath
ifeq ($(orthufuserversion),2.2.6)
	@echo "orthofuser right version is already installed"
else
	@echo "version ${orthufuserversion}"
	@echo "orthofuser is installed, but not the right version"
	cd ${DIR}/software/OrthoFinder/ && git pull
endif
else
	@echo "orthofuser is not installed and needs to be installed"
	cd ${DIR}/software && git clone https://github.com/macmanes-lab/OrthoFinder.git
	@echo PATH=\$$PATH:${DIR}/software/OrthoFinder/orthofinder >> pathfile
endif

postscript:brew setup shmlast lastal mcl samtools hmmer quorum orthofuser rcorrector blast spades trinity shannon seqtk busco trimmomatic transrate bowtie2 salmon
	@printf "\n\n*** The following location(s), if any print, need to be added to your PATH ***"
	@printf "\n*** They will be automatically to your ~/.profile or ~/.bash_profile ***\n\n"
	@cat pathfile
	@cat pathfile >> ~/.profile
	@cat pathfile >> ~/.bash_profile
	@export PATH=$$PATH:$$(cat pathfile)
	@printf "\n\n\n"
	@printf "\n*** type <<source ~/.profile>> to complete the install ***\n\n"
