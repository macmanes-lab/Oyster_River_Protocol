#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	oyster.mk prep main READ1= READ2= CPU=24 RUNOUT=runname
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
RUNOUT =
ASSEMBLY=
LINEAGE=
BUSCOUT := BUSCO_$(shell basename ${ASSEMBLY} .fasta)
BUSCODB :=
START=1
INPUT := $(shell basename ${READ1})


prep: setup run_scripts
main: run_trimmomatic run_rcorrector run_trinity run_spades run_shannon orthofusing report
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

run_trimmomatic:
	trimmomatic PE -threads $(CPU) -baseout ${DIR}/rcorr/${RUNOUT}.TRIM.fastq ${READ1} ${READ2}  LEADING:3 TRAILING:3 ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:30:10 MINLEN:25


run_rcorrector:
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 31 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.fastq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.fastq -od ${DIR}/rcorr


run_trinity:
	cd ${DIR}/assemblies && \
	Trinity --no_normalize_reads --seqType fq --output ${RUNOUT}.trinity --max_memory 50G --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup

run_spades:
	cd ${DIR}/assemblies && \
	rnaspades.py --only-assembler -o ${RUNOUT}.spades_k75 --threads $(CPU) --memory 100 -k 75 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq && \
	rnaspades.py --only-assembler -o ${RUNOUT}.spades_k55 --threads $(CPU) --memory 100 -k 55 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq && \
	mv ${RUNOUT}.spades_k55/transcripts.fasta ${RUNOUT}.transcripts55.fasta && \
	mv ${RUNOUT}.spades_k75/transcripts.fasta ${RUNOUT}.transcripts75.fasta  && \
	rm -fr ${RUNOUT}.spades_k55 ${RUNOUT}.spades_k75

run_shannon:
	cd ${DIR}/assemblies && \
	python $$(which shannon.py) -o ${RUNOUT}.shannon --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -p $(CPU) -K 75 && \
	mv ${RUNOUT}.shannon/shannon.fasta ${RUNOUT}.shannon.fasta && \
	rm -fr ${RUNOUT}.shannon

orthofusing:
	cd ${DIR}/orthofuse && \
	mkdir -p ${RUNOUT} && \
	ln -s ${DIR}/assemblies/${RUNOUT}.transcripts55.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.transcripts55.fasta && \
	ln -s ${DIR}/assemblies/${RUNOUT}.transcripts75.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.transcripts75.fasta && \
	ln -s ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.trinity.Trinity.fasta && \
	ln -s ${DIR}/assemblies/${RUNOUT}.shannon.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.shannon.fasta && \
	python $$(which orthofuser.py) -f ${DIR}/orthofuse/${RUNOUT}/ -og -t $(CPU) -a $(CPU) && \
	cat ${DIR}/orthofuse/${RUNOUT}/*fasta > ${DIR}/orthofuse/${RUNOUT}/merged.fasta && \
	transrate -o ${DIR}/orthofuse/${RUNOUT}/merged -t $(CPU) -a ${DIR}/orthofuse/${RUNOUT}/merged.fasta --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq && \
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${RUNOUT}/ -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}') && \
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${RUNOUT}/ -name Orthogroups.txt 2> /dev/null) && \
	for i in $$(eval echo "{1..$$END}") ; do sed -n ''$$i'p' $$ORTHOINPUT | tr ' ' '\n' > ${DIR}/orthofuse/${RUNOUT}/$$i.groups; done && \
	echo All the text files are made, start GREP  && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *groups 2> /dev/null | parallel -j $(CPU) "grep -wf {} $$(find ${DIR}/orthofuse/${RUNOUT}/ -name contigs.csv 2> /dev/null) > {1}.orthout 2> /dev/null" && \
	echo About to delete all the text files  && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *groups -delete && \
	echo Search output files  && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *orthout 2> /dev/null | parallel -j $(CPU) "awk -F, 'BEGIN {max = 0} {if (\$$9>max) max=\$$9} END {print \$$1 \"\\t\" max}'" | tee -a ${DIR}/orthofuse/${RUNOUT}/good.list && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *orthout -delete && \
	python $$(which filter.py) ${DIR}/orthofuse/${RUNOUT}/merged.fasta <(awk '{print $$1}' ${DIR}/orthofuse/${RUNOUT}/good.list) > ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta && \
	touch orthofuse.done

busco.done:
	cd ${DIR}/reports && \
	python $$(which run_BUSCO.py) -i ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta -m transcriptome --cpu $(CPU) -l /mnt/lustre/macmaneslab/macmanes/BUSCODB/${LINEAGE} -o ${BUSCOUT} && \
	touch busco.done

transrate.done:
	cd ${DIR}/reports && \
	#transrate -o transrate_${basename ${RUNOUT}/orthomerged.fasta .fasta}  -a ${DIR}/orthofuse/${RUNOUT}/orthomerged.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate -o transrate_${basename ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta .fasta}  -a ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -t $(CPU) && \
	#transrate -o transrate_${basename ${ASSEMBLY} .fasta}  -a ${DIR}/assemblies/${ASSEMBLY} --left ${DIR}/rcorr/${RUNOUT}.skewer-trimmed-pair1.fastq --right ${DIR}/rcorr/${RUNOUT}.skewer-trimmed-pair2.fastq -t $(CPU) && \
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
