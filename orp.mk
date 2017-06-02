#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	for i in 1 2 5 10 20 40 60 80 100; do ./orp.mk prep main SAMP=$i CPU=24; done
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
RCORR ?= ${shell which rcorrector}
RCORRDIR := $(dir $(firstword $(RCORR)))
READ1=
READ2=
BUSCO ?= ${shell which run_BUSCO.py}
BUSCODIR := $(dir $(firstword $(BUSCO)))
DATASET := $(shell basename ${READ1} _1.fastq.gz)
ASSEMBLY=
LINEAGE=
BUSCOUT := BUSCO_$(shell basename ${ASSEMBLY} .fasta)
BUSCODB :=


prep: setup run_scripts
main: subsamp_reads run_rcorrector run_skewer rcorr_trinity rcorr_spades rcorr_shannon transfuse
report:busco.done transrate.done report
busco:busco.done
transrate:transrate.done

.DELETE_ON_ERROR:
.PHONY:report transfuse

setup:
	mkdir -p ${DIR}/scripts
	mkdir -p ${DIR}/reads
	mkdir -p ${DIR}/assemblies
	mkdir -p ${DIR}/rcorr
	mkdir -p ${DIR}/reports

run_scripts:
	@echo Downloading Scripts
	cd ${DIR}/scripts && \
	curl -LO https://raw.githubusercontent.com/macmanes-lab/general/master/filter.py && \
	wget https://raw.githubusercontent.com/macmanes/read_error_corr/master/barcodes.fa

subsamp_reads:
	cd ${DIR}/reads && \
	seqtk sample -s102340 ${READ1} ${SAMP}000000 > ${DATASET}.${SAMP}.subsamp_1.fastq && \
	seqtk sample -s102340 ${READ2} ${SAMP}000000 > ${DATASET}.${SAMP}.subsamp_2.fastq

run_rcorrector:
	cd ${DIR}/rcorr && \
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 31 -1 ${DIR}/reads/${DATASET}.${SAMP}.subsamp_1.fastq -2 ${DIR}/reads/${DATASET}.${SAMP}.subsamp_2.fastq && \
	awk -F 'l:' '{print $$1}' ${DIR}/rcorr/${DATASET}.${SAMP}.subsamp_1.cor.fq | sed 's_ __g' > tmp && mv tmp ${DIR}/rcorr/${DATASET}.${SAMP}.subsamp_1.cor.fq && \
	awk -F 'l:' '{print $$1}' ${DIR}/rcorr/${DATASET}.${SAMP}.subsamp_2.cor.fq | sed 's_ __g' > tmp && mv tmp ${DIR}/rcorr/${DATASET}.${SAMP}.subsamp_2.cor.fq

run_skewer:
	cd ${DIR}/rcorr && \
	skewer -l 25 -m pe -o ${DATASET}.${SAMP}.skewer --mean-quality 2 --end-quality 2 -t $(CPU) -x ${DIR}/scripts/barcodes.fa ${DIR}/rcorr/${DATASET}.${SAMP}.subsamp_1.cor.fq ${DIR}/rcorr/${DATASET}.${SAMP}.subsamp_2.cor.fq

rcorr_trinity:
	cd ${DIR}/assemblies && \
	Trinity --no_normalize_reads --seqType fq --output ${DATASET}.${SAMP}.trinity --max_memory 50G --left ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup

rcorr_spades:
	cd ${DIR}/assemblies && \
	rnaspades.py -o ${DATASET}.${SAMP}.spades_k75 --threads $(CPU) --memory 100 -k 75 -1 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq -2 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq && \
	rnaspades.py -o ${DATASET}.${SAMP}.spades_k55 --threads $(CPU) --memory 100 -k 55 -1 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq -2 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq && \
	mv ${DATASET}.${SAMP}.spades_k55/transcripts.fasta ${DATASET}.${SAMP}.transcripts55.fasta && \
	mv ${DATASET}.${SAMP}.spades_k75/transcripts.fasta ${DATASET}.${SAMP}.transcripts75.fasta  && \
	rm -fr ${DATASET}.${SAMP}.spades_k55 ${DATASET}.${SAMP}.spades_k75

rcorr_shannon:
	cd ${DIR}/assemblies && \
	python $$(which shannon.py) -o ${DATASET}.${SAMP}.shannon --left ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq -p $(CPU) -K 75 && \
	mv ${DATASET}.${SAMP}.shannon/shannon.fasta ${DATASET}.${SAMP}.shannon.fasta && \
	rm -fr ${DATASET}.${SAMP}.shannon

transfuse:
	cd ${DIR}/assemblies && \
	transfuse -t $(CPU) -i 0.98 -o ${DATASET}.${SAMP}.transfuse.fasta -l ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq -r ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq -a ${DATASET}.${SAMP}.transcripts55.fasta,${DATASET}.${SAMP}.transcripts75.fasta,${DATASET}.${SAMP}.trinity.Trinity.fasta,${DATASET}.${SAMP}.shannon.fasta

busco.done:
	cd ${DIR}/reports && \
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${ASSEMBLY} -m transcriptome --cpu $(CPU) -l /mnt/lustre/macmaneslab/macmanes/BUSCODB/${LINEAGE} -o ${BUSCOUT} && \
	touch busco.done

transrate.done:
	cd ${DIR}/reports && \
	transrate -o transrate_${basename ${ASSEMBLY} .fasta}  -a ${DIR}/assemblies/${ASSEMBLY} --left ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq -t $(CPU) && \
	touch transrate.done

report:
	printf "\n\n*****  QUALITY REPORT FOR: ${ASSEMBLY} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/run_${BUSCOUT} -name short*) | sed -n 8p  | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
	printf " \n\n"
