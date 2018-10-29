#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#       quant.mk all MEM= CPU= JOBS= SAMPLE= SUFFIX= TRANSCRIPTOME=
#

###############################################################
## More detailed information about the job
###############################################################

VERSION = 2.0.0 # suggest we move to 2.1.0
MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
MEM=128
JOBS=1
RCORR := ${shell which rcorrector}
RCORRDIR := $(dir $(firstword $(RCORR)))
TRANSCRIPTOME=
SAMPLE=
SUFFIX=
START=1
INPUT := $(shell basename ${READ1})  ## NOT NEEDED?
brewpath := $(shell which brew 2>/dev/null)
rcorrpath := $(shell which rcorrector 2>/dev/null)
trimmomaticpath := $(shell which trimmomatic 2>/dev/null)
salmonpath := $(shell which salmon 2>/dev/null)


check:
setup:${DIR}/quant_setup.done
salmon:${DIR}/${SAMPLE}salmonquant.done

all: check setup salmon

.DELETE_ON_ERROR:
.PHONY:report check

${DIR}/quant_setup.done:
	@mkdir -p ${DIR}/reports
	@mkdir -p ${DIR}/quants
	touch ${DIR}/quant_setup.done

check:
ifdef salmonpath
else
    $(error "\n\n*** SALMON is not installed, must fix ***")
endif
ifdef trimmomaticpath
else
    $(error  "\n\n Maybe TRIMMOMATIC is not installed, or maybe you are working on Bridges")
endif
ifdef rcorrpath
else
    $(error "\n\n *** RCORRECTOR is not installed, must fix ***")
endif

${DIR}/${SAMPLE}salmonquant.done:
	salmon index --no-version-check -t ${TRANSCRIPTOME}  -i ${TRANSCRIPTOME}.ortho.idx --type quasi -k 31
	salmon quant --no-version-check -p ${CPU} -i ${TRANSCRIPTOME}.ortho.idx --seqBias --gcBias -l a -1 ${DIR}/${SAMPLE}1${SUFFIX} -2 ${DIR}/${SAMPLE}2${SUFFIX} -o ${DIR}/quants/${SAMPLE};\
	touch ${DIR}/${SAMPLE}salmonquant.done

	printf " \n\n"
	echo ${SAMPLE} quantification complete
	runtime=$(shell time)
	echo $(runtime)
	source ${MAKEDIR}/software/anaconda/install/bin/deactivate 
