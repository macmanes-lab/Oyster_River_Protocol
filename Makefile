#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	make
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
# The conda this installer bootstraps when there is none on PATH. A conda
# already on PATH (a system install, a cluster module, Miniforge) is used as
# it is instead, and never updated or reconfigured.
BUNDLED_CONDA := ${DIR}/software/anaconda/install
ORP_ENVS := orp orp_spades orp_trinity orp_busco orp_transabyss
SPROT_URL := https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
diamond_data := $(and $(wildcard ${DIR}/software/diamond/swissprot.dmnd),$(wildcard ${DIR}/software/diamond/uniprot_sprot.fasta))
busco_data := $(shell find ${DIR}/busco_dbs -iname "eukaryota_odb12*" -type d 2>/dev/null)
VERSION := ${shell cat  ${MAKEDIR}version.txt}
# Everything Trinity needs to build and run, but not Trinity itself: the
# bioconda package ships a ParaFly built without OpenMP (see
# scripts/build_trinity.sh). cmake<4 because Inchworm and Chrysalis declare a
# cmake_minimum_required that cmake 4 refuses. bwa, seqtk and bashplotlib are for strandeval.
TRINITY_ENV_PACKAGES = compilers "cmake<4" make git kmer-jellyfish bowtie2 samtools salmon=1.10.3 "openjdk>=17" perl perl-db_file python numpy bwa=0.7.19 bashplotlib seqtk=1.5 libgomp zlib
TRINITY_DIR := ${DIR}/software/trinityrnaseq
trinity_built := $(shell nm -D ${TRINITY_DIR}/trinity-plugins/BIN/ParaFly 2>/dev/null | grep -c GOMP_ | grep -v '^0$$')

# Start of every recipe that uses conda: find it at run time, not when this
# file is parsed, so a conda the `conda` target has just bootstrapped is
# found too. Sets conda_exe, and solver to mamba when the base has one (it is
# faster to solve with) or else conda, whose libmamba solver is the default
# from 23.10.
WITH_CONDA = set -e; \
	conda_exe=$$(command -v conda || echo "${BUNDLED_CONDA}/bin/conda"); \
	conda_base=$$("$$conda_exe" info --base); \
	if [ -x "$$conda_base/bin/mamba" ]; then solver="$$conda_base/bin/mamba"; \
	elif command -v mamba >/dev/null; then solver=mamba; \
	else solver="$$conda_exe"; fi

.PHONY: all setup conda orp trinity diamond_data busco_data postscript clean

all: setup conda orp trinity diamond_data busco_data postscript

.DELETE_ON_ERROR:

setup:
	@mkdir -p ${DIR}/shared
	@mkdir -p ${DIR}/software/anaconda
	@mkdir -p ${DIR}/software/diamond

conda:setup
	@if command -v conda >/dev/null; then \
		echo "conda: using $$(command -v conda)"; \
	elif [ -x ${BUNDLED_CONDA}/bin/conda ]; then \
		echo "conda: using ${BUNDLED_CONDA}/bin/conda"; \
	else \
		set -e; \
		cd ${DIR}/software/anaconda; \
		curl -fLO https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh; \
		bash Anaconda3-2020.11-Linux-x86_64.sh -b -p ${BUNDLED_CONDA}; \
		${BUNDLED_CONDA}/bin/conda update -y -n base conda; \
		${BUNDLED_CONDA}/bin/conda install -y -n base -c conda-forge mamba; \
		echo ". ${BUNDLED_CONDA}/etc/profile.d/conda.sh" > ${DIR}/pathfile; \
		echo "PATH=\$$PATH:${BUNDLED_CONDA}/bin" >> ${DIR}/pathfile; \
	fi

# One env at a time, each created only if no env has exactly its name, so a
# partial install picks up where it stopped.
orp:orp_env.yml conda setup
	@$(WITH_CONDA); \
	envs=$$("$$conda_exe" env list | awk '!/^#/ {print $$1}'); \
	created=0; \
	have_env() { echo "$$envs" | grep -x "$$1" >/dev/null || { created=1; false; }; }; \
	have_env orp_spades || "$$solver" create -y -c bioconda -c conda-forge --override-channels --name orp_spades spades=4.3.0 python=3.14; \
	have_env orp_trinity || "$$solver" create -y -c bioconda -c conda-forge --override-channels --name orp_trinity ${TRINITY_ENV_PACKAGES}; \
	have_env orp_busco || "$$solver" create -y -c bioconda -c conda-forge --override-channels --name orp_busco busco=6.1.0; \
	have_env orp_transabyss || "$$solver" create -y -c bioconda -c conda-forge --override-channels --name orp_transabyss transabyss=2.0.1; \
	have_env orp || "$$solver" env create -f ${DIR}/orp_env.yml; \
	if [ "$$created" = 1 ] && [ "$$conda_base" = "${BUNDLED_CONDA}" ]; then "$$solver" clean -ya; fi; \
	echo "conda envs: ${ORP_ENVS} present"


trinity:orp
ifdef trinity_built
		@echo "trinity is already built"
else
		@$(WITH_CONDA); \
		"$$conda_exe" run --no-capture-output -n orp_trinity ${DIR}/scripts/build_trinity.sh ${TRINITY_DIR}
endif

# Downloaded and unzipped via .part files, so an interrupted download is
# fetched again rather than taken for finished.
diamond_data:orp
ifdef diamond_data
		@echo "diamond_data is already installed"
else
		@$(WITH_CONDA); \
		cd ${DIR}/software/diamond; \
		if [ ! -s uniprot_sprot.fasta ]; then \
			curl -fL -o uniprot_sprot.fasta.gz.part ${SPROT_URL}; \
			gzip -dc uniprot_sprot.fasta.gz.part > uniprot_sprot.fasta.part; \
			mv uniprot_sprot.fasta.part uniprot_sprot.fasta; \
			rm -f uniprot_sprot.fasta.gz.part; \
		fi; \
		"$$conda_exe" run -n orp diamond makedb --in uniprot_sprot.fasta -d swissprot
endif

busco_data:orp
ifdef busco_data
		@echo "busco_data is already installed"
else
		@$(WITH_CONDA); \
		mkdir -p ${DIR}/busco_dbs; \
		"$$conda_exe" run -n orp_busco busco --download eukaryota_odb12.2 --download_path ${DIR}/busco_dbs
endif

# pathfile exists only when the conda target bootstrapped conda. Each of its
# lines goes into ~/.profile -- and ~/.bash_profile if there is one: creating
# it would stop bash reading ~/.profile -- unless the file already has it, so
# re-running make adds nothing twice. ~/.bashrc gets the conda.sh line, for
# shells that are not login shells.
postscript: setup orp trinity diamond_data busco_data conda
	@if [ -f pathfile ]; then\
		printf "\n\n*** The following location(s), if any print, need to be added to your PATH ***";\
		printf "\n*** They will be automatically to your ~/.profile or ~/.bash_profile ***\n\n";\
		cat pathfile;\
		while IFS= read -r line; do\
			for rc in ~/.profile ~/.bash_profile; do\
				if [ "$$rc" = ~/.bash_profile ] && [ ! -f "$$rc" ]; then continue; fi;\
				grep -qxF "$$line" "$$rc" 2>/dev/null || echo "$$line" >> "$$rc";\
			done;\
		done < pathfile;\
		conda_sh=$$(head -n 1 pathfile);\
		grep -qxF "$$conda_sh" ~/.bashrc 2>/dev/null || echo "$$conda_sh" >> ~/.bashrc;\
		printf "\n\n\n";\
		printf "\n*** type ``source ~/.profile`` to complete the install ***\n\n";\
	fi

# Removes the five ORP envs from whichever conda holds them, then the
# bootstrapped conda (if any) and Trinity's build.
clean:
	-@conda_exe=$$(command -v conda || echo "${BUNDLED_CONDA}/bin/conda"); \
	for env in ${ORP_ENVS}; do "$$conda_exe" remove -y -n $$env --all 2>/dev/null || true; done
	rm -fr ${DIR}/software/anaconda/install
	rm -fr ${DIR}/software/transabyss
	rm -fr ${TRINITY_DIR}
	rm -fr ${DIR}/software/anaconda/
	rm -fr ${DIR}/pathfile
