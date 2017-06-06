#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	oyster.mk prep main READ1= READ2= CPU=24
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
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
main: run_rcorrector run_skewer rcorr_trinity rcorr_spades rcorr_shannon orthofusing report
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
	cd ${DIR}/scripts && \
	curl -LO https://raw.githubusercontent.com/macmanes-lab/general/master/filter.py && \
	wget https://raw.githubusercontent.com/macmanes/read_error_corr/master/barcodes.fa

run_rcorrector:
	cd ${DIR}/rcorr && \
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 31 -1 ../${READ1} -2 ../${READ1} && \
	awk -F 'l:' '{print $$1}' ${DIR}/rcorr/${DATASET}.1.cor.fq | sed 's_ __g' > tmp && mv tmp ${DIR}/rcorr/${DATASET}.1.cor.fq && \
	awk -F 'l:' '{print $$1}' ${DIR}/rcorr/${DATASET}.2.cor.fq | sed 's_ __g' > tmp && mv tmp ${DIR}/rcorr/${DATASET}.2.cor.fq

run_skewer:
	cd ${DIR}/rcorr && \
	skewer -l 25 -m pe -o ${DATASET}.skewer --mean-quality 2 --end-quality 2 -t $(CPU) -x ${DIR}/scripts/barcodes.fa ${DIR}/rcorr/${DATASET}.1.cor.fq ${DIR}/rcorr/${DATASET}.2.cor.fq

rcorr_trinity:
	cd ${DIR}/assemblies && \
	Trinity --no_normalize_reads --seqType fq --output ${DATASET}.trinity --max_memory 50G --left ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair2.fastq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup

rcorr_spades:
	cd ${DIR}/assemblies && \
	rnaspades.py -o ${DATASET}.spades_k75 --threads $(CPU) --memory 100 -k 75 -1 ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair1.fastq -2 ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair2.fastq && \
	rnaspades.py -o ${DATASET}.spades_k55 --threads $(CPU) --memory 100 -k 55 -1 ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair1.fastq -2 ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair2.fastq && \
	mv ${DATASET}.spades_k55/transcripts.fasta ${DATASET}.transcripts55.fasta && \
	mv ${DATASET}.spades_k75/transcripts.fasta ${DATASET}.transcripts75.fasta  && \
	rm -fr ${DATASET}.spades_k55 ${DATASET}.spades_k75

rcorr_shannon:
	cd ${DIR}/assemblies && \
	python $$(which shannon.py) -o ${DATASET}.shannon --left ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair2.fastq -p $(CPU) -K 75 && \
	mv ${DATASET}.shannon/shannon.fasta ${DATASET}.shannon.fasta && \
	rm -fr ${DATASET}.shannon

orthofusing:
	cd ${DIR}/orthofuse && \
	mkdir ${DATASET} && \
	ln -s ${DIR}/assemblies/${DATASET}.transcripts55.fasta ${DIR}/orthofuse/${DATASET}/${DATASET}.transcripts55.fasta && \
	ln -s ${DIR}/assemblies/${DATASET}.transcripts75.fasta ${DIR}/orthofuse/${DATASET}/${DATASET}.transcripts75.fasta && \
	ln -s ${DIR}/assemblies/${DATASET}.trinity.Trinity.fasta ${DIR}/orthofuse/${DATASET}/${DATASET}.trinity.Trinity.fasta && \
	ln -s ${DIR}/assemblies/${DATASET}.shannon.fasta ${DIR}/orthofuse/${DATASET}/${DATASET}.shannon.fasta && \
	python $$(which orthofinder.py) -f ${DIR}/orthofuse/${DATASET}/ -og -t $(CPU) -a $(CPU) && \
	cat ${DIR}/orthofuse/${DATASET}/*fasta > ${DIR}/orthofuse/${DATASET}/merged.fasta && \
	transrate -o ${DIR}/orthofuse/${DATASET}/merged -t $(CPU) -a ${DIR}/orthofuse/${DATASET}/merged.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} && \
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${DATASET}/ -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}') && \
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${DATASET}/ -name Orthogroups.txt 2> /dev/null) && \
	for i in $$(eval echo "{1..$$END}") ; do sed -n ''$$i'p' $$ORTHOINPUT | tr ' ' '\n' | grep -f - $$(find ${DIR}/orthofuse/${DATASET}/ -name contigs.csv 2> /dev/null) | awk -F, 'BEGIN {max = 0} {if ($$9>max) max=$$9} END {print $$1 "\t" max}' | tee -a ${DIR}/orthofuse/${DATASET}/good.list; done && \
	python $$(which filter.py) ${DIR}/orthofuse/${DATASET}/merged.fasta <(awk '{print $$1}' ${DATASET}/good.list) > ${DIR}/orthofuse/${DATASET}/${DATASET}.orthomerged.fasta && \
	touch orthofuse.done

busco.done:
	cd ${DIR}/reports && \
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${ASSEMBLY} -m transcriptome --cpu $(CPU) -l /mnt/lustre/macmaneslab/macmanes/BUSCODB/${LINEAGE} -o ${BUSCOUT} && \
	touch busco.done

transrate.done:
	cd ${DIR}/reports && \
	#transrate -o transrate_${basename ${DATASET}/orthomerged.fasta .fasta}  -a ${DIR}/orthofuse/${DATASET}/orthomerged.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate -o transrate_${basename ${ASSEMBLY} .fasta}  -a ${DIR}/assemblies/${ASSEMBLY} --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	#transrate -o transrate_${basename ${ASSEMBLY} .fasta}  -a ${DIR}/assemblies/${ASSEMBLY} --left ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${DATASET}.skewer-trimmed-pair2.fastq -t $(CPU) && \
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
