#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	oyster.mk prep main READ1= READ2= MEM=500 CPU=24 RUNOUT=runname
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
MEM=128
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

run_trimmomatic:
run_rcorrector:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
run_trinity:${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta
run_spades55:${DIR}/assemblies/${RUNOUT}.spades55.fasta
run_spades75:${DIR}/assemblies/${RUNOUT}.spades75.fasta
run_shannon:${DIR}/assemblies/${RUNOUT}.shannon.fasta
merge:${DIR}/orthofuse/${RUNOUT}/merged.fasta
orthofusing:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta

main: setup run_trimmomatic run_rcorrector run_trinity run_spades75 run_spades55 run_shannon merge orthofusing report
report:busco.done transrate.done reportgen
busco:busco.done
transrate:transrate.done

.DELETE_ON_ERROR:
.PHONY:report

setup:
	@mkdir -p ${DIR}/scripts
	@mkdir -p ${DIR}/reads
	@mkdir -p ${DIR}/assemblies
	@mkdir -p ${DIR}/rcorr
	@mkdir -p ${DIR}/reports
	@mkdir -p ${DIR}/orthofuse

${DIR}/rcorr/${RUNOUT}.TRIM_1P.fastq:
	trimmomatic PE -threads $(CPU) -baseout ${DIR}/rcorr/${RUNOUT}.TRIM.fastq ${READ1} ${READ2} LEADING:3 TRAILING:3 ILLUMINACLIP:${MAKEDIR}/barcodes/barcodes.fa:2:30:10 MINLEN:25

${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq:${DIR}/rcorr/${RUNOUT}.TRIM_1P.fastq
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 31 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.fastq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.fastq -od ${DIR}/rcorr

${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	Trinity --no_normalize_reads --seqType fq --output ${DIR}/assemblies/${RUNOUT}.trinity --max_memory $(MEM)G --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup

${DIR}/assemblies/${RUNOUT}.spades55.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	rnaspades.py --only-assembler -o ${DIR}/assemblies/${RUNOUT}.spades_k55 --threads $(CPU) --memory $(MEM) -k 55 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq
	mv ${DIR}/assemblies/${RUNOUT}.spades_k55/transcripts.fasta ${DIR}/assemblies/${RUNOUT}.spades55.fasta
	rm -fr ${DIR}/assemblies/${RUNOUT}.spades_k55

${DIR}/assemblies/${RUNOUT}.spades75.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	rnaspades.py --only-assembler -o ${DIR}/assemblies/${RUNOUT}.spades_k75 --threads $(CPU) --memory $(MEM) -k 75 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq
	mv ${DIR}/assemblies/${RUNOUT}.spades_k75/transcripts.fasta ${DIR}/assemblies/${RUNOUT}.spades75.fasta
	rm -fr ${DIR}/assemblies/${RUNOUT}.spades_k75

${DIR}/assemblies/${RUNOUT}.shannon.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	seqtk seq -A ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq > /mouse/orthofinder/orp/rcorr/$$(basename ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq .fq).fa
	seqtk seq -A ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq > /mouse/orthofinder/orp/rcorr/$$(basename ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq .fq).fa
	python $$(which shannon.py) -o ${DIR}/assemblies/${RUNOUT}.shannon --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fa --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fa -p $(CPU) -K 75
	mv ${DIR}/assemblies/${RUNOUT}.shannon/shannon.fasta {DIR}/assemblies/${RUNOUT}.shannon.fasta
	rm -fr {DIR}/assemblies/${RUNOUT}.shannon

${DIR}/orthofuse/${RUNOUT}/merged.fasta:${DIR}/assemblies/${RUNOUT}.spades55.fasta ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta ${DIR}/assemblies/${RUNOUT}.shannon.fasta
	cd ${DIR}/orthofuse && \
	mkdir -p ${DIR}/orthofuse/${RUNOUT} && \
	ln -sf ${DIR}/assemblies/${RUNOUT}.spades55.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.spades55.fasta && \
	ln -sf ${DIR}/assemblies/${RUNOUT}.spades75.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.spades75.fasta && \
	ln -sf ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.trinity.Trinity.fasta && \
	ln -sf ${DIR}/assemblies/${RUNOUT}.shannon.fasta ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.shannon.fasta && \
	python $$(which orthofuser.py) -I 4 -f ${DIR}/orthofuse/${RUNOUT}/ -og -t $(CPU) -a $(CPU) && \
	cat ${DIR}/orthofuse/${RUNOUT}/*fasta > ${DIR}/orthofuse/${RUNOUT}/merged.fasta

${DIR}/assemblies/${RUNOUT}.orthomerged.fasta:${DIR}/orthofuse/${RUNOUT}/merged.fasta
	cd ${DIR}/orthofuse && \
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${RUNOUT}/ -cmin 1 -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}') && \
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${RUNOUT}/ -cmin 1 -name Orthogroups.txt 2> /dev/null) && \
	parallel  -j $(CPU) -k "sed -n ''$$i'p' $$ORTHOINPUT | tr ' ' '\n' > ${DIR}/orthofuse/${RUNOUT}/{1}.groups"  ::: $$(eval echo "{1..$$END}") && \
	transrate -o ${DIR}/orthofuse/${RUNOUT}/merged -t $(CPU) -a ${DIR}/orthofuse/${RUNOUT}/merged.fasta --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq && \
	echo All the text files are made, start GREP  && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *groups 2> /dev/null | parallel -j $(CPU) "grep -wf {} $$(find ${DIR}/orthofuse/${RUNOUT}/ -name contigs.csv 2> /dev/null) > {1}.orthout 2> /dev/null" && \
	echo About to delete all the text files  && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *groups -delete && \
	echo Search output files  && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *orthout 2> /dev/null | parallel -j $(CPU) "awk -F, 'BEGIN {max = 0} {if (\$$9>max) max=\$$9} END {print \$$1 \"\\t\" max}'" | tee -a ${DIR}/orthofuse/${RUNOUT}/good.list && \
	find ${DIR}/orthofuse/${RUNOUT}/ -name *orthout -delete && \
	python $$(which filter.py) ${DIR}/orthofuse/${RUNOUT}/merged.fasta <(awk '{print $$1}' ${DIR}/orthofuse/${RUNOUT}/good.list) > ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta && \
	cp ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta

busco.done:
	cd ${DIR}/reports && \
	python $$(which run_BUSCO.py) -i ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta -m transcriptome --cpu $(CPU) -l /mnt/lustre/macmaneslab/macmanes/BUSCODB/${LINEAGE} -o ${BUSCOUT} && \
	touch busco.done

transrate.done:
	cd ${DIR}/reports && \
	transrate -o transrate_${basename ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta .fasta}  -a ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -t $(CPU) && \
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
