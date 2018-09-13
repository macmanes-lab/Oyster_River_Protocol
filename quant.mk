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
# RUNOUT =
TRANSCRIPTOME=
# READS=
SAMPLE=
SUFFIX=
START=1
INPUT := $(shell basename ${READ1})  ## NOT NEEDED?
# FASTADIR=
brewpath := $(shell which brew 2>/dev/null)
rcorrpath := $(shell which rcorrector 2>/dev/null)
trimmomaticpath := $(shell which trimmomatic 2>/dev/null)
salmonpath := $(shell which salmon 2>/dev/null)
seqtkpath := $(shell which seqtk 2>/dev/null)    ### needed?
#THREADS := $(${CPUS} / ${JOBS})    ###### add the math here (cpus/jobs)


check:
setup:${DIR}/quant_setup.done
run_trimmomatic:${DIR}/${SAMPLE}trimmomatic.done
run_rcorrector:${DIR}/${SAMPLE}rcorr.done
salmon:${DIR}/${SAMPLE}salmonquant.done

all: check setup run_trimmomatic run_rcorrector salmon

.DELETE_ON_ERROR:
.PHONY:report check

${DIR}/quant_setup.done:
	@mkdir -p ${DIR}/rcorr
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

${DIR}/${SAMPLE}trimmomatic.done:
	@if [ $$(hostname | cut -d. -f3-5) == 'bridges.psc.edu' ];\
	then\
		java -jar $$TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -threads ${CPU} -baseout ${DIR}/rcorr/${SAMPLE}TRIM.fastq ${DIR}/${SAMPLE}1${SUFFIX} ${DIR}/${SAMPLE}2${SUFFIX} LEADING:3 TRAILING:3 ILLUMINACLIP:${MAKEDIR}/barcodes/barcodes.fa:2:30:10 MINLEN:25;\
		touch ${DIR}/${SAMPLE}trimmomatic.done;\
	else\
		trimmomatic PE -threads ${CPU} -baseout ${DIR}/rcorr/${SAMPLE}TRIM.fastq ${DIR}/${SAMPLE}1${SUFFIX} ${DIR}/${SAMPLE}2${SUFFIX} LEADING:3 TRAILING:3 ILLUMINACLIP:${MAKEDIR}/barcodes/barcodes.fa:2:30:10 MINLEN:25;\
		touch ${DIR}/${SAMPLE}trimmomatic.done;\
	fi

${DIR}/${SAMPLE}rcorr.done:${DIR}/${SAMPLE}trimmomatic.done
	perl ${RCORRDIR}/run_rcorrector.pl -t ${CPU} -k 31 -1 ${DIR}/rcorr/${SAMPLE}TRIM_1P.fastq -2 ${DIR}/rcorr/${SAMPLE}TRIM_2P.fastq -od ${DIR}/rcorr;\
	touch ${DIR}/${SAMPLE}rcorr.done

${DIR}/${SAMPLE}salmonquant.done:${DIR}/${SAMPLE}rcorr.done
	salmon index --no-version-check -t ${DIR}/${TRANSCRIPTOME}  -i ${TRANSCRIPTOME}.ortho.idx --type quasi -k 31
	salmon quant --no-version-check -p ${CPU} -i ${TRANSCRIPTOME}.ortho.idx --seqBias --gcBias -l a -1 ${DIR}/rcorr/${SAMPLE}TRIM_1P.cor.fq -2 ${DIR}/rcorr/${SAMPLE}TRIM_2P.cor.fq -o ${DIR}/quants/{};\
	rm -fr ${RUNOUT}.ortho.idx
	touch ${DIR}/${SAMPLE}salmonquant.done

	printf " \n\n"
	echo ${SAMPLE} quantification complete
	source ${MAKEDIR}/software/anaconda/install/bin/deactivate orp_v2
