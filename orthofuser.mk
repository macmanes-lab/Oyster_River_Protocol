#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	orthofuser.mk all READ1= READ2= CPU= RUNOUT= FASTADIR= LINEAGE=
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
READ1=
READ2=
BUSCO := ${shell which run_BUSCO.py}
BUSCODIR := $(dir $(firstword $(BUSCO)))
RUNOUT =
LINEAGE=
BUSCODB :=
INPUT := $(shell basename ${READ1})
FASTADIR=


merge:${DIR}/orthofuse/${RUNOUT}/merged.fasta
orthotransrate:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
orthofusing:${FASTADIR}/${RUNOUT}.orthomerged.fasta

all: setup merge orthotransrate orthofusing report
busco:${DIR}/reports/busco.done
transrate:${DIR}/reports/transrate.done
report:busco transrate reportgen



.DELETE_ON_ERROR:
.PHONY:report

setup:
	@mkdir -p ${DIR}/reports

${DIR}/orthofuse/${RUNOUT}/merged.fasta:
	mkdir -p ${DIR}/orthofuse/${RUNOUT}/working
	for fasta in $$(ls ${FASTADIR}); do python ${MAKEDIR}/scripts/long.seq.py ${FASTADIR}/$$fasta ${DIR}/orthofuse/${RUNOUT}/working/$$fasta.short.fasta 200; done
	python $$(which orthofuser.py) -I 4 -f ${DIR}/orthofuse/${RUNOUT}/working/ -og -t $(CPU) -a $(CPU)
	cat ${DIR}/orthofuse/${RUNOUT}/working/*short.fasta > ${DIR}/orthofuse/${RUNOUT}/merged.fasta

${DIR}/orthofuse/${RUNOUT}/orthotransrate.done:${DIR}/orthofuse/${RUNOUT}/merged.fasta
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${RUNOUT}/working/ -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}') && \
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${RUNOUT}/working/ -name Orthogroups.txt 2> /dev/null) && \
	echo $$(eval echo "{1..$$END}") | tr ' ' '\n' > list && \
	cat list | parallel  -j $(CPU) -k "sed -n ''{}'p' $$ORTHOINPUT | tr ' ' '\n' | sed '1d' > ${DIR}/orthofuse/${RUNOUT}/{1}.groups"
	transrate -o ${DIR}/orthofuse/${RUNOUT}/merged -t $(CPU) -a ${DIR}/orthofuse/${RUNOUT}/merged.fasta --left ${READ1} --right ${READ2}
	touch ${DIR}/orthofuse/${RUNOUT}/orthotransrate.done

${FASTADIR}/${RUNOUT}.orthomerged.fasta:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
	echo All the text files are made, start GREP
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*groups' 2> /dev/null | parallel -j $(CPU) "grep -wf {} $$(find ${DIR}/orthofuse/${RUNOUT}/ -name contigs.csv 2> /dev/null) > {1}.orthout 2> /dev/null"
	echo About to delete all the text files
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*groups' -delete
	echo Search output files
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' 2> /dev/null | parallel -j $(CPU) "awk -F, -v max=0 '{if(\$$14>max){want=\$$1; max=\$$14}}END{print want}'" >> ${DIR}/orthofuse/${RUNOUT}/good.list
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' -delete
	python ${MAKEDIR}/scripts/filter.py ${DIR}/orthofuse/${RUNOUT}/merged.fasta ${DIR}/orthofuse/${RUNOUT}/good.list > ${FASTADIR}/${RUNOUT}.orthomerged.fasta
	rm ${DIR}/orthofuse/${RUNOUT}/good.list

${DIR}/reports/busco.done:${FASTADIR}/${RUNOUT}.orthomerged.fasta
	python $$(which run_BUSCO.py) -i ${FASTADIR}/${RUNOUT}.orthomerged.fasta -m transcriptome --cpu $(CPU) -o ${RUNOUT}
	mv run_${RUNOUT} ${DIR}/reports/
	touch ${DIR}/reports/busco.done

${DIR}/reports/transrate.done:${FASTADIR}/${RUNOUT}.orthomerged.fasta
	transrate -o ${DIR}/reports/transrate_${RUNOUT}  -a ${FASTADIR}/${RUNOUT}.orthomerged.fasta --left ${READ1} --right ${READ2} -t $(CPU)
	touch ${DIR}/reports/transrate.done

${DIR}/quants/orthomerged/quant.sf:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta  -i ${RUNOUT}.ortho.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.ortho.idx --seqBias --gcBias -l a -1 ${READ1} -2 ${READ1} -o ${DIR}/quants/salmon_orthomerged_${RUNOUT}
	rm -fr ${RUNOUT}.ortho.idx

reportgen:
	printf "\n\n*****  QUALITY REPORT FOR: ${RUNOUT} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/run_${RUNOUT} -name 'short*') | sed -n 8p  | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf " \n\n"
