#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	make
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CONDAROOT = ${DIR}/software/anaconda/install/
transrate := $(shell ls ${DIR}/software/orp-transrate/transrate 2>/dev/null)
diamond_data := $(shell ls ${DIR}/software/diamond/uniprot_sprot.fasta 2>/dev/null)
busco_data := $(shell find ${DIR}/busco_dbs -iname "eukaryota_odb12*" -type d 2>/dev/null)
conda := $(shell conda info 2>/dev/null)
orp := $(shell ${DIR}/software/anaconda/install/bin/conda info --envs | grep orp 2>/dev/null)
VERSION := ${shell cat  ${MAKEDIR}version.txt}

all: setup conda orp transrate diamond_data busco_data postscript

.DELETE_ON_ERROR:

setup:
	@mkdir -p ${DIR}/shared
	@mkdir -p ${DIR}/software/anaconda
	@mkdir -p ${DIR}/software/diamond

conda:setup
ifdef conda
else
	cd ${DIR}/software/anaconda && curl -LO https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
	cd ${DIR}/software/anaconda && bash Anaconda3-2020.11-Linux-x86_64.sh -b -p install/
	@echo ". ${DIR}/software/anaconda/install/etc/profile.d/conda.sh" >> ~/.bashrc;
	@echo ". ${DIR}/software/anaconda/install/etc/profile.d/conda.sh" > pathfile;
	source ~/.bashrc;
endif

orp:orp_env.yml conda setup
ifdef orp
else
	( \
				source ${DIR}/software/anaconda/install/etc/profile.d/conda.sh; \
				conda activate; \
				conda update -y -n base conda; \
				conda config --add channels conda-forge; \
				conda config --add channels bioconda; \
				conda install mamba -n base -yc conda-forge; \
				mamba create -y -c bioconda -c conda-forge --override-channels --name orp_spades spades=4.3.0 python=3.14; \
				mamba create -y -c bioconda -c conda-forge --override-channels --name orp_trinity trinity=2.15.2 bwa=0.7.19 bashplotlib seqtk=1.5 salmon=1.10.3; \
				mamba create -y -c bioconda -c conda-forge --override-channels --name orp_busco busco=6.1.0; \
				mamba create -y -c bioconda -c conda-forge --override-channels --name orp_transabyss transabyss=2.0.1; \
				mamba create -y -c bioconda -c conda-forge --override-channels --name orp_orthofinder orthofinder=3.1.5; \
				mamba env create -f ${DIR}/orp_env.yml; \
				mamba clean -ya; \
				conda deactivate; \
  )
	@echo PATH=\$$PATH:${DIR}/software/anaconda/install/bin >> pathfile;
endif


diamond_data:conda
ifdef diamond_data
		@echo "diamond_data is already installed"
else
		cd ${DIR}/software/diamond && curl -LO ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && gzip -d uniprot_sprot.fasta.gz && ${DIR}/software/anaconda/install/envs/orp/bin/diamond makedb --in uniprot_sprot.fasta -d swissprot
endif

busco_data:conda
ifdef busco_data
		@echo "busco_data is already installed"
else
		mkdir -p ${DIR}/busco_dbs
		${DIR}/software/anaconda/install/envs/orp_busco/bin/busco --download eukaryota_odb12.2 --download_path ${DIR}/busco_dbs
endif

transrate:
ifdef transrate
else
	cd ${DIR}/software && tar -zxf orp-transrate.tar.gz
	@echo PATH=\$$PATH:${DIR}/software/orp-transrate >> pathfile
endif

postscript: setup orp diamond_data busco_data conda transrate
	@if [ -f pathfile ]; then\
		printf "\n\n*** The following location(s), if any print, need to be added to your PATH ***";\
		printf "\n*** They will be automatically to your ~/.profile or ~/.bash_profile ***\n\n";\
		cat pathfile;\
		cat pathfile >> ~/.profile;\
		cat pathfile >> ~/.bash_profile;\
		cat pathfile >> ~/.bash_profile;\
		export PATH=$$PATH:$$(cat pathfile);\
		printf "\n\n\n";\
		printf "\n*** type ``source ~/.profile`` to complete the install ***\n\n";\
	fi

clean:
	${DIR}/software/anaconda/install/bin/conda remove -y --name orp --all
	rm -fr ${DIR}/software/anaconda/install
	rm -fr ${DIR}/software/orp-transrate
	rm -fr ${DIR}/software/transabyss
	rm -fr ${DIR}/software/anaconda/
	rm -fr ${DIR}/pathfile
