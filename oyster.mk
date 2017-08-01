#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	oyster.mk main READ1= READ2= MEM=500 CPU=24 RUNOUT=runname
# oyster.mk orthofuse FASTADIR= READ1= READ2= MEM=500 CPU=24 RUNOUT=runname
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
FASTADIR=

run_trimmomatic:
run_rcorrector:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
run_trinity:${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta
run_spades55:${DIR}/assemblies/${RUNOUT}.spades55.fasta
run_spades75:${DIR}/assemblies/${RUNOUT}.spades75.fasta
run_shannon:${DIR}/assemblies/${RUNOUT}.shannon.fasta
merge:${DIR}/orthofuse/${RUNOUT}/merged.fasta
orthotransrate:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
orthofusing:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta


main: setup run_trimmomatic run_rcorrector run_trinity run_spades75 run_spades55 run_shannon merge orthotransrate orthofusing report
report:busco transrate reportgen
busco:${DIR}/reports/busco.done
transrate:${DIR}/reports/transrate.done

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
	@if [ $$(hostname | cut -d. -f3-5) == 'bridges.psc.edu' ];\
	then\
		java -jar $$TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -threads $(CPU) -baseout ${DIR}/rcorr/${RUNOUT}.TRIM.fastq ${READ1} ${READ2} LEADING:3 TRAILING:3 ILLUMINACLIP:${MAKEDIR}/barcodes/barcodes.fa:2:30:10 MINLEN:25;\
	else\
		trimmomatic PE -threads $(CPU) -baseout ${DIR}/rcorr/${RUNOUT}.TRIM.fastq ${READ1} ${READ2} LEADING:3 TRAILING:3 ILLUMINACLIP:${MAKEDIR}/barcodes/barcodes.fa:2:30:10 MINLEN:25;\
	fi

${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq:${DIR}/rcorr/${RUNOUT}.TRIM_1P.fastq
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 31 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.fastq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.fastq -od ${DIR}/rcorr

${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	Trinity --no_version_check --no_normalize_reads --seqType fq --output ${DIR}/assemblies/${RUNOUT}.trinity --max_memory $(MEM)G --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup
	awk '{print $$1}' ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta > ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fa && mv -f ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fa ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta

${DIR}/assemblies/${RUNOUT}.spades55.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	rnaspades.py --only-assembler -o ${DIR}/assemblies/${RUNOUT}.spades_k55 --threads $(CPU) --memory $(MEM) -k 55 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq
	mv ${DIR}/assemblies/${RUNOUT}.spades_k55/transcripts.fasta ${DIR}/assemblies/${RUNOUT}.spades55.fasta
	rm -fr ${DIR}/assemblies/${RUNOUT}.spades_k55

${DIR}/assemblies/${RUNOUT}.spades75.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	rnaspades.py --only-assembler -o ${DIR}/assemblies/${RUNOUT}.spades_k75 --threads $(CPU) --memory $(MEM) -k 75 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq
	mv ${DIR}/assemblies/${RUNOUT}.spades_k75/transcripts.fasta ${DIR}/assemblies/${RUNOUT}.spades75.fasta
	rm -fr ${DIR}/assemblies/${RUNOUT}.spades_k75

${DIR}/assemblies/${RUNOUT}.shannon.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	seqtk seq -A ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq > ${DIR}/rcorr/$$(basename ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq .fq).fa
	seqtk seq -A ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq > ${DIR}/rcorr/$$(basename ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq .fq).fa
	python $$(which shannon.py) -o ${DIR}/assemblies/${RUNOUT}.shannon --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fa --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fa -p $(CPU) -K 75
	mv ${DIR}/assemblies/${RUNOUT}.shannon/shannon.fasta ${DIR}/assemblies/${RUNOUT}.shannon.fasta
	rm -fr {DIR}/assemblies/${RUNOUT}.shannon

${DIR}/orthofuse/${RUNOUT}/merged.fasta:
	mkdir -p ${DIR}/orthofuse/${RUNOUT}/working
	for fasta in $$(ls ${DIR}/assemblies/${RUNOUT}*); do python ${MAKEDIR}/scripts/long.seq.py ${DIR}/assemblies/$$fasta ${DIR}/orthofuse/${RUNOUT}/working/$$fasta.short.fasta 200; done
	python $$(which orthofuser.py) -I 4 -f ${DIR}/orthofuse/${RUNOUT}/working/ -og -t $(CPU) -a $(CPU)
	cat ${DIR}/orthofuse/${RUNOUT}/working/*short.fasta > ${DIR}/orthofuse/${RUNOUT}/merged.fasta

${DIR}/orthofuse/${RUNOUT}/orthotransrate.done:${DIR}/orthofuse/${RUNOUT}/merged.fasta
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${RUNOUT}/working/ -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}') && \
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${RUNOUT}/working/ -name Orthogroups.txt 2> /dev/null) && \
	echo $$(eval echo "{1..$$END}") | tr ' ' '\n' > list && \
	cat list | parallel  -j $(CPU) -k "sed -n ''{}'p' $$ORTHOINPUT | tr ' ' '\n' | sed '1d' > ${DIR}/orthofuse/${RUNOUT}/{1}.groups"
	transrate -o ${DIR}/orthofuse/${RUNOUT}/merged -t $(CPU) -a ${DIR}/orthofuse/${RUNOUT}/merged.fasta --left ${READ1} --right ${READ2}
	touch ${DIR}/orthofuse/${RUNOUT}/orthotransrate.done

${DIR}/assemblies/${RUNOUT}.orthomerged.fasta:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
	echo All the text files are made, start GREP
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*groups' 2> /dev/null | parallel -j $(CPU) "grep -wf {} $$(find ${DIR}/orthofuse/${RUNOUT}/ -name contigs.csv 2> /dev/null) > {1}.orthout 2> /dev/null"
	echo About to delete all the text files
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*groups' -delete
	echo Search output files
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' 2> /dev/null | parallel -j $(CPU) "awk -F, -v max=0 '{if(\$$14>max){want=\$$1; max=\$$14}}END{print want}'" | tee -a ${DIR}/orthofuse/${RUNOUT}/good.list
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' -delete
	python ${MAKEDIR}/scripts/filter.py ${DIR}/orthofuse/${RUNOUT}/merged.fasta ${DIR}/orthofuse/${RUNOUT}/good.list > ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta
	cp ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	rm ${DIR}/orthofuse/${RUNOUT}/good.list

${DIR}/reports/busco.done:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	python3 $$(which run_BUSCO.py) -i ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta -m transcriptome --cpu $(CPU) -l ${LINEAGE} -o ${RUNOUT}
	mv run_${RUNOUT} ${DIR}/reports/
	touch ${DIR}/reports/busco.done

${DIR}/reports/transrate.done:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	transrate -o ${DIR}/reports/transrate_${RUNOUT}  -a ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta --left ${READ1} --right ${READ2} -t $(CPU)
	touch ${DIR}/reports/transrate.done

reportgen:
	printf "\n\n*****  QUALITY REPORT FOR: ${RUNOUT} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/run_${RUNOUT} -name 'short*') | sed -n 8p  | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf " \n\n"
