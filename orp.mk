#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	for i in 1 2 5 10 20 40 60 80 100; do ./correct.mk main SAMP=$i CPU=36; done
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
RCORR ?= ${shell which rcorrector}
RCORRDIR := $(dir $(firstword $(RCORR)))
READ1=
READ2=
L = $(basename ${READ1} .fastq)
R = $(basename ${READ2} .fastq)
BUSCO ?= ${shell which BUSCO_v1.22.py}
BUSCODIR := $(dir $(firstword $(BUSCO)))
ASSEMBLY=
LINEAGE=
BUSCOUT := BUSCO_$(basename ${ASSEMBLY} .fasta)



all: prep main
prep: setup run_scripts
main: run_rcorrector run_skewer rcorr_trinity rcorr_binpacker transfuse
report:busco.done transrate.done report
busco:busco.done
transrate:transrate.done

.DELETE_ON_ERROR:
.PHONY:report

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

run_rcorrector:
	cd ${DIR}/rcorr && \
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 31 -1 ${DIR}/reads/${READ1} -2 ${DIR}/reads/${READ2}

run_skewer:
	cd ${DIR}/rcorr && \
	skewer -l 25 -m pe -o skewer --mean-quality 2 --end-quality 2 -t $(CPU) -x ${DIR}/scripts/barcodes.fa ${DIR}/rcorr/$(join ${L},.cor.fq) ${DIR}/rcorr/$(join ${R},.cor.fq)

rcorr_trinity:
	cd ${DIR}/assemblies && \
	Trinity --seqType fq --output trinity_rcorr31 --max_memory 50G --left ${DIR}/rcorr/skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/skewer-trimmed-pair2.fastq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup

rcorr_binpacker:
	cd ${DIR}/assemblies && \
	BinPacker -d -q -s fq -p pair -m RF -k 25 -g 200 -l ${DIR}/rcorr/skewer-trimmed-pair1.fastq -r ${DIR}/rcorr/skewer-trimmed-pair2.fastq

transfuse:
	cd ${DIR}/assemblies && \
	transfuse -t $(CPU) -i 0.98 -o transfuse -l ${DIR}/rcorr/skewer-trimmed-pair1.fastq -r ${DIR}/rcorr/skewer-trimmed-pair2.fastq -a Rcorr_binpacker/BinPacker.fa,trinity_rcorr31.Trinity.fasta

busco.done:
    cd ${DIR}/reports && \
    python3 ${BUSCODIR}BUSCO_v1.22.py -in ${DIR}/assemblies/${ASSEMBLY} -m trans --cpu $(CPU) -l ${BUSCODIR}${LINEAGE} -o ${BUSCOUT} && \
    touch busco.done

transrate.done:
    cd ${DIR}/reports && \
    transrate -o transrate_${basename ${ASSEMBLY} .fasta}  -a ${DIR}/assemblies/${ASSEMBLY} --left ${DIR}/rcorr/skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/skewer-trimmed-pair2.fastq -t $(CPU) && \
    touch transrate.done

report:
    printf "\n\n*****  QUALITY REPORT FOR: ${ASSEMBLY} **** \n\n"
    printf "*****  BUSCO SCORE ~~~~~>           " | tee qualreport.${basename ${ASSEMBLY} .fasta}
    cat $$(find reports/run_${BUSCOUT} -name short*) | sed -n 5p  | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
    printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
    cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$41}' | sed -n 2p | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
    printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
    cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$42}' | sed -n 2p | tee -a qualreport.${basename ${ASSEMBLY} .fasta}
    printf " \n\n"
