#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	oyster.mk main READ1= READ2= MEM=110 CPU=24 RUNOUT=runname STRAND=RF
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
BUSCO_THREADS=${CPU}
MEM=110
READ1=
READ2=
ASSEMBLY=
RUNOUT =
START=1
STRAND :=
TPM_FILT =
FASTADIR=
salmonpath := $(shell which salmon 2>/dev/null)
VERSION := ${shell cat  ${MAKEDIR}version.txt}

.DEFAULT_GOAL := main

help:
main: setup check welcome readcheck cdhit salmon filter
preprocess:setup check welcome readcheck run_trimmomatic run_rcorrector
update_merge:setup check welcome readcheck merge orthotransrate orthofusing diamond posthack cdhit salmon filter busco transrate strandeval report
salmon:${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf
diamond:${DIR}/assemblies/diamond/${RUNOUT}.trinity.diamond.txt
clean:
setup:${DIR}/assemblies/working ${DIR}/reads ${DIR}/rcorr ${DIR}/assemblies/diamond ${DIR}/assemblies ${DIR}/reports ${DIR}/orthofuse ${DIR}/quants
cdhit:${DIR}/assemblies/${RUNOUT}.ORP.fasta
filter:${DIR}/assemblies/${RUNOUT}.filter.done

.DELETE_ON_ERROR:
.PHONY:report check clean

${DIR}/assemblies/working ${DIR}/reads ${DIR}/rcorr ${DIR}/assemblies/diamond ${DIR}/assemblies ${DIR}/reports ${DIR}/orthofuse ${DIR}/quants:
	@mkdir -p ${DIR}/reads
	@mkdir -p ${DIR}/assemblies
	@mkdir -p ${DIR}/reports
	@mkdir -p ${DIR}/quants
	@mkdir -p ${DIR}/assemblies/diamond
	@mkdir -p ${DIR}/assemblies/working

check:
ifdef salmonpath
else
	$(error "\n\n*** SALMON is not installed, must fix ***")
endif



help:
	printf "\n\n*****  Welcome to the Oyster River Prptocol ***** \n"
	printf "*****  This is version ${VERSION} *****\n\n"
	printf "Usage:\n\n"
	printf "/path/to/Oyster_River/Protocol/oyster.mk main CPU=24 \n"
	printf "MEM=128\n"
	printf "STRAND=RF\n"
	printf "READ1=1.subsamp_1.cor.fq\n"
	printf "READ2=1.subsamp_2.cor.fq\n"
	printf "RUNOUT=test\n\n"


welcome:
	printf "\n\n*****  Welcome to the Oyster River ***** \n"
	printf "*****  This is the filtering script, version ${VERSION} ***** \n\n "
	printf " \n\n"

readcheck:
	if [ -e ${READ1} ]; then printf ""; else printf "\n\n\n\n ERROR: YOUR READ1 FILE DOES NOT EXIST AT THE LOCATION YOU SPECIFIED\n\n\n\n "; $$(shell exit); fi;
	if [ -e ${READ2} ]; then printf ""; else printf "\n\n\n\n ERROR: YOUR READ2 FILE DOES NOT EXIST AT THE LOCATION YOU SPECIFIED\n\n\n\n "; $$(shell exit); fi;
ifeq ($(shell file ${READ1} | awk '{print $$2}'),gzip)
	if [ $$(gzip -cd $${READ1} | head -n 400 | awk '{if(NR%4==2) {count++; bases += length} } END{print int(bases/count)}') -gt 75 ] && [ $$(gzip -cd $${READ2} | head -n 400 | awk '{if(NR%4==2) {count++; bases += length} } END{print int(bases/count)}') -gt 75 ];\
	then\
		printf " ";\
	else\
		printf "\n\n\n\n IT LOOKS LIKE YOUR READS ARE NOT AT LEAST 76 BP LONG,\n ";\
		printf "PLEASE EDIT YOUR COMMAND USING THE "SPADES2_KMER=INT" FLAGS,\n";\
		printf " SETTING THE ASSEMBLY KMER LENGTH TO AN ODD NUMBER LESS THAN YOUR READ LENGTH \n\n\n\n";\
		$$(shell exit);\
	fi
else
	if [ $$(head -n400 $${READ1} | awk '{if(NR%4==2) {count++; bases += length} } END{print int(bases/count)}') -gt 75 ] && [ $$(head -n400 $${READ2} | awk '{if(NR%4==2) {count++; bases += length} } END{print int(bases/count)}') -gt 75 ];\
	then\
		printf " ";\
	else\
		printf "\n\n\n\n IT LOOKS LIKE YOUR READS ARE NOT AT LEAST 75 BP LONG,\n ";\
		printf "PLEASE EDIT YOUR COMMAND USING THE "SPADES2_KMER=INT" FLAGS,\n";\
		printf " SETTING THE ASSEMBLY KMER LENGTH LESS THAN YOUR READ LENGTH \n\n\n\n";\
		$$(shell exit);\
	fi
endif



${DIR}/assemblies/${RUNOUT}.ORP.fasta:${ASSEMBLY}
	cd ${DIR}/assemblies/ && cd-hit-est -M 5000 -T $(CPU) -c .98 -i ${ASSEMBLY} -o ${DIR}/assemblies/${RUNOUT}.ORP.fasta
	diamond blastx -p $(CPU) -e 1e-8 --top 0.1 -q ${DIR}/assemblies/${RUNOUT}.ORP.fasta -d ${MAKEDIR}/software/diamond/swissprot  -o ${DIR}/assemblies/${RUNOUT}.ORP.diamond.txt
	awk '{print $$2}' ${DIR}/assemblies/${RUNOUT}.ORP.diamond.txt | awk -F "|" '{print $$3}' | cut -d _ -f2 | sort | uniq | wc -l > ${DIR}/assemblies/working/${RUNOUT}.unique.ORP.txt
	rm ${DIR}/assemblies/${RUNOUT}.ORP.fasta.clstr

${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf:${DIR}/assemblies/${RUNOUT}.ORP.fasta ${READ1} ${READ2}
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.ORP.fasta  -i ${RUNOUT}.ortho.idx --type quasi -k 31
	salmon quant --no-version-check --validateMappings -p $(CPU) -i ${RUNOUT}.ortho.idx --seqBias --gcBias -l a -1 ${READ1} -2 ${READ1} -o ${DIR}/quants/salmon_orthomerged_${RUNOUT}
	rm -fr ${RUNOUT}.ortho.idx

${DIR}/assemblies/${RUNOUT}.filter.done ${DIR}/assemblies/working/${RUNOUT}.saveme.fasta:${DIR}/assemblies/${RUNOUT}.ORP.fasta ${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf
ifdef TPM_FILT
	cat ${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf| awk '$$4 > $(TPM_FILT)' | cut -f1 | sed 1d > ${DIR}/assemblies/working/${RUNOUT}.HIGHEXP.txt
	cat ${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf| awk '$$4 < $(TPM_FILT)' | cut -f1 | sed 1d > ${DIR}/assemblies/working/${RUNOUT}.LOWEXP.txt
	cp ${DIR}/assemblies/${RUNOUT}.ORP.fasta ${DIR}/assemblies/working/${RUNOUT}.ORP_BEFORE_TPM_FILT.fasta
	python ${MAKEDIR}/scripts/filter.py ${DIR}/assemblies/${RUNOUT}.ORP.fasta ${DIR}/assemblies/working/${RUNOUT}.HIGHEXP.txt > ${DIR}/assemblies/working/${RUNOUT}.ORP.HIGHEXP.fasta
	grep -Fwf ${DIR}/assemblies/working/${RUNOUT}.LOWEXP.txt ${DIR}/assemblies/${RUNOUT}.ORP.diamond.txt >> ${DIR}/assemblies/working/${RUNOUT}.blasted
	awk '{print $$1}' ${DIR}/assemblies/working/${RUNOUT}.blasted | sort | uniq | tee -a ${DIR}/assemblies/working/${RUNOUT}.donotremove.list
	python ${MAKEDIR}/scripts/filter.py ${DIR}/assemblies/${RUNOUT}.ORP.fasta ${DIR}/assemblies/working/${RUNOUT}.donotremove.list > ${DIR}/assemblies/working/${RUNOUT}.saveme.fasta
	cat ${DIR}/assemblies/working/${RUNOUT}.saveme.fasta ${DIR}/assemblies/working/${RUNOUT}.ORP.HIGHEXP.fasta > ${DIR}/assemblies/working/${RUNOUT}.tmp.fasta
	mv ${DIR}/assemblies/working/${RUNOUT}.tmp.fasta ${DIR}/assemblies/${RUNOUT}.ORP.fasta
	touch ${DIR}/assemblies/${RUNOUT}.filter.done
else
	touch ${DIR}/assemblies/${RUNOUT}.filter.done
endif



clean:
	rm -fr ${DIR}/orthofuse/${RUNOUT}/ ${DIR}/rcorr/${RUNOUT}.TRIM_2P.fastq ${DIR}/rcorr/${RUNOUT}.TRIM_1P.fastq ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta \
	${DIR}/assemblies/${RUNOUT}.spades55.fasta ${DIR}/assemblies/${RUNOUT}.spades75.fasta ${DIR}/assemblies/${RUNOUT}.transabyss.fasta \
	${DIR}/orthofuse/${RUNOUT}/merged.fasta ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta ${DIR}/assemblies/diamond/${RUNOUT}.trinity.diamond.txt \
	${DIR}/assemblies/diamond/${RUNOUT}.newbies.fasta ${DIR}/reports/busco.done ${DIR}/reports/transrate.done ${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf \
	${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq ${DIR}/reports/run_${RUNOUT}.orthomerged/ ${DIR}/reports/transrate_${RUNOUT}/
