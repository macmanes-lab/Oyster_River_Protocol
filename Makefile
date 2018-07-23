#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	make
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CONDAROOT = ${DIR}/software/anaconda/install/
orthopath := $(shell source ${DIR}/software/anaconda/install/bin/activate py27; which orthofuser.py 2>/dev/null; conda deactivate)
orthufuserversion = $(shell source ${DIR}/software/anaconda/install/bin/activate py27; orthofuser.py --help | grep "OrthoFinder version" | awk '{print $$3}';conda deactivate)
transrate := $(shell which transrate 2>/dev/null)
transabysspath := $(shell which ${DIR}/software/transabyss/transabyss 2>/dev/null)
transabyssversion = $(shell transabyss --version)
diamond_data := $(shell ls ${DIR}/software/diamond/uniprot_sprot.fasta 2>/dev/null)
busco_data := $(shell ls ${DIR}/busco_dbs/eukaryota_odb9 2>/dev/null)
conda := $(shell conda info 2>/dev/null)
orp_v2 := $(shell conda info --envs | grep orp_v2 2>/dev/null)
py27 := $(shell conda info --envs | grep py27 2>/dev/null)


all: setup conda environment orthofuser transrate transabyss diamond_data busco_data postscript

.DELETE_ON_ERROR:

setup:
	@mkdir -p ${DIR}/shared
	@mkdir -p ${DIR}/software/anaconda
	@mkdir -p ${DIR}/software/diamond
	@rm -f pathfile

conda:
ifdef conda
else
	cd ${DIR}/software/anaconda && curl -LO https://repo.anaconda.com/archive/Anaconda3-5.1.0-Linux-x86_64.sh
	cd ${DIR}/software/anaconda && bash Anaconda3-5.1.0-Linux-x86_64.sh -b -p ${DIR}/software/anaconda/install
endif

orp_v2:environment.yml conda
ifdef orp_v2
else
	( \
       source ${DIR}/software/anaconda/install/bin/activate; \
       conda update -y -n base conda; \
			 source ${DIR}/software/anaconda/install/bin/deactivate; \
			 conda env create -f environment.yml python=3.6; \
  )
	@echo PATH=\$$PATH:${DIR}/software/anaconda/install/bin >> pathfile;
endif

py27:conda
ifdef py27
else
	( \
       source ${DIR}/software/anaconda/install/bin/activate; \
       conda update -y -n base conda; \
			 source ${DIR}/software/anaconda/install/bin/deactivate; \
			 conda create -y -n py27 python=2.7 anaconda; \
  )
	@echo PATH=\$$PATH:${DIR}/software/anaconda/install/bin >> pathfile;
endif

transabyss:
ifdef transabysspath
ifeq ($(transabyssversion),2.0.1)
	@echo "TransABySS is already installed"
else
	@echo "version ${transabyssversion}"
	@echo "TransABySS is installed, but not the right version"
	cd ${DIR}/software/transabyss && git pull
endif
else
	cd ${DIR}/software/ && git clone https://github.com/bcgsc/transabyss.git
	@echo PATH=\$$PATH:${DIR}/software/transabyss >> pathfile
endif

diamond_data:conda
ifdef diamond_data
	@echo "diamond_data is already installed"
else
	 cd ${DIR}/software/diamond && curl -LO ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && gzip -d uniprot_sprot.fasta.gz
	 cd ${DIR}/software/diamond && ${DIR}/software/anaconda/install/envs/orp_v2/bin/diamond makedb --in uniprot_sprot.fasta -d swissprot
endif

busco_data:
ifdef busco_data
else
	mkdir ${DIR}/busco_dbs && cd ${DIR}/busco_dbs
	cd ${DIR}/busco_dbs && wget http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz && tar -zxf eukaryota_odb9.tar.gz
endif

transrate:
	cd ${DIR}/software && tar -zxf orp-transrate.tar.gz
	@echo PATH=\$$PATH:${DIR}/software/orp-transrate >> pathfile

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

postscript: setup diamond_data busco_data orthofuser conda transrate
	@printf "\n\n*** The following location(s), if any print, need to be added to your PATH ***"
	@printf "\n*** They will be automatically to your ~/.profile or ~/.bash_profile ***\n\n"
	@cat pathfile
	@cat pathfile >> ~/.profile
	@cat pathfile >> ~/.bash_profile
	@export PATH=$$PATH:$$(cat pathfile)
	source $$HOME/.profile
	@printf "\n\n\n"
	@printf "\n*** type <<source ~/.profile>> to complete the install ***\n\n"
