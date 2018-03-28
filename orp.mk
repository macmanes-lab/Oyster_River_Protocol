#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	for i in 1 2 5 10 20 40 60 80 100; do ./orp.mk prep main SAMP=$i CPU=24; done
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
MEM=120
RCORR := ${shell which rcorrector}
RCORRDIR := $(dir $(firstword $(RCORR)))
READ1=
READ2=
BUSCO := ${shell which run_BUSCO.py}
BUSCODIR := $(dir $(firstword $(BUSCO)))
DATASET := $(shell basename ${READ1} _1.fastq.gz)
ASSEMBLY=
LINEAGE=
BUSCOUT := BUSCO_$(shell basename ${ASSEMBLY} .fasta)
BUSCODB :=
START=1

prep: setup run_scripts
main: subsamp_reads run_rcorrector run_skewer rcorr_trinity rcorr_spades rcorr_shannon orthofusing
report:busco.done transrate.done reportgen
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
	mkdir -p ${DIR}/orthofuse

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
	rnaspades.py -o ${DATASET}.${SAMP}.spades_k75 --memory $(MEM) --threads $(CPU) -k 75 -1 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq -2 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq && \
	rnaspades.py -o ${DATASET}.${SAMP}.spades_k55 --memory $(MEM) --threads $(CPU) -k 55 -1 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq -2 ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq && \
	mv ${DATASET}.${SAMP}.spades_k55/transcripts.fasta ${DATASET}.${SAMP}.transcripts55.fasta && \
	mv ${DATASET}.${SAMP}.spades_k75/transcripts.fasta ${DATASET}.${SAMP}.transcripts75.fasta  && \
	rm -fr ${DATASET}.${SAMP}.spades_k55 ${DATASET}.${SAMP}.spades_k75

rcorr_shannon:
	cd ${DIR}/assemblies && \
	python $$(which shannon.py) -o ${DATASET}.${SAMP}.shannon --left ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${DATASET}.${SAMP}.skewer-trimmed-pair2.fastq -p $(CPU) -K 75 && \
	mv ${DATASET}.${SAMP}.shannon/shannon.fasta ${DATASET}.${SAMP}.shannon.fasta && \
	rm -fr ${DATASET}.${SAMP}.shannon

orthofusing:
	cd ${DIR}/orthofuse && \
	mkdir -p ${DATASET}.${SAMP} && \
	ln -s ${DIR}/assemblies/${DATASET}.${SAMP}.transcripts55.fasta ${DIR}/orthofuse/${DATASET}.${SAMP}/${DATASET}.${SAMP}.transcripts55.fasta && \
	ln -s ${DIR}/assemblies/${DATASET}.${SAMP}.transcripts75.fasta ${DIR}/orthofuse/${DATASET}.${SAMP}/${DATASET}.${SAMP}.transcripts75.fasta && \
	ln -s ${DIR}/assemblies/${DATASET}.${SAMP}.trinity.Trinity.fasta ${DIR}/orthofuse/${DATASET}.${SAMP}/${DATASET}.${SAMP}.trinity.Trinity.fasta && \
	ln -s ${DIR}/assemblies/${DATASET}.${SAMP}.shannon.fasta ${DIR}/orthofuse/${DATASET}.${SAMP}/${DATASET}.${SAMP}.shannon.fasta && \
	python $$(which orthofinder.py) -f ${DIR}/orthofuse/${DATASET}.${SAMP}/ -og -t $(CPU) -a $(CPU) && \
	cat ${DIR}/orthofuse/${DATASET}.${SAMP}/*fasta > ${DIR}/orthofuse/${DATASET}.${SAMP}/merged.fasta && \
	transrate -o ${DIR}/orthofuse/${DATASET}.${SAMP}/merged -t $(CPU) -a ${DIR}/orthofuse/${DATASET}.${SAMP}/merged.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} && \
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${DATASET}.${SAMP}/ -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}') && \
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${DATASET}.${SAMP}/ -name Orthogroups.txt 2> /dev/null) && \
	for i in $$(eval echo "{1..$$END}") ; do sed -n ''$$i'p' $$ORTHOINPUT | tr ' ' '\n' > ${DIR}/orthofuse/${DATASET}.${SAMP}/$$i.groups; done && \
	echo All the text files are made, start GREP  && \
	find ${DIR}/orthofuse/${DATASET}.${SAMP}/ -name *groups 2> /dev/null | parallel -j $(CPU) "grep -wf {} $$(find ${DIR}/orthofuse/${DATASET}.${SAMP}/ -name contigs.csv 2> /dev/null) > {1}.orthout 2> /dev/null" && \
	echo About to delete all the text files  && \
	find ${DIR}/orthofuse/${DATASET}.${SAMP}/ -name *groups -delete && \
	echo Search output files  && \
	find ${DIR}/orthofuse/${DATASET}.${SAMP}/ -name *orthout 2> /dev/null | parallel -j $(CPU) "awk -F, 'BEGIN {max = 0} {if (\$$9>max) max=\$$9} END {print \$$1 \"\\t\" max}'" | tee -a ${DIR}/orthofuse/${DATASET}.${SAMP}/good.list && \
	find ${DIR}/orthofuse/${DATASET}.${SAMP}/ -name *orthout -delete && \
	python $$(which filter.py) ${DIR}/orthofuse/${DATASET}.${SAMP}/merged.fasta <(awk '{print $$1}' ${DIR}/orthofuse/${DATASET}.${SAMP}/good.list) > ${DIR}/orthofuse/${DATASET}.${SAMP}/${DATASET}.${SAMP}.orthomerged.fasta && \
	touch orthofuse.done

busco.done:
	cd ${DIR}/reports && \
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${ASSEMBLY} --force -m transcriptome --cpu $(CPU) -l /mnt/lustre/macmaneslab/macmanes/BUSCODB/${LINEAGE} -o ${BUSCOUT} && \
	touch busco.done

transrate.done:
	cd ${DIR}/reports && \
	transrate -o transrate_$(basename ${ASSEMBLY} .fasta)  -a ${DIR}/assemblies/${ASSEMBLY} --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	touch transrate.done

reportgen:
	printf "\n\n*****  QUALITY REPORT FOR: ${ASSEMBLY} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/run_${BUSCOUT} -name short*) | sed -n 8p  | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	printf " \n\n"
