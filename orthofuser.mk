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
orthofusing:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta

all: setup merge orthotransrate orthofusing report
report:busco.done transrate.done reportgen
busco:busco.done
transrate:transrate.done


.DELETE_ON_ERROR:
.PHONY:report

setup:
	@mkdir -p ${DIR}/reports
	@mkdir -p ${DIR}/orthofuse

${DIR}/orthofuse/${RUNOUT}/merged.fasta:
	mkdir -p ${DIR}/orthofuse/${RUNOUT}
	cd ${FASTADIR} && \
	for fasta in $$(ls $(basename ${FASTADIR})); do python ${MAKEDIR}/scripts/long.seq.py $$fasta $$($$fasta).short.fasta;gzip $$fasta; done
	python $$(which orthofuser.py) -I 4 -f ${FASTADIR} -og -t $(CPU) -a $(CPU)
	mv ${FASTADIR}/Results* ${DIR}/orthofuse/${RUNOUT}/
	cat ${FASTADIR}/* > ${DIR}/orthofuse/${RUNOUT}/merged.fasta

${DIR}/orthofuse/${RUNOUT}/orthotransrate.done:${DIR}/orthofuse/${RUNOUT}/merged.fasta
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${RUNOUT}/ -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}')
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${RUNOUT}/ -name Orthogroups.txt 2> /dev/null)
	parallel  -j $(CPU) -k "sed -n ''{}'p' $$ORTHOINPUT | tr ' ' '\n' | sed '1d' > ${DIR}/orthofuse/${RUNOUT}/{1}.groups"  ::: $$(eval echo "{1..$$END}")
	transrate -o ${DIR}/orthofuse/${RUNOUT}/merged -t $(CPU) -a ${DIR}/orthofuse/${RUNOUT}/merged.fasta --left ${READ1} --right ${READ2}
	touch ${DIR}/orthofuse/${RUNOUT}/orthotransrate.done

${DIR}/assemblies/${RUNOUT}.orthomerged.fasta:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
	echo All the text files are made, start GREP
	find ${DIR}/orthofuse/${RUNOUT}/ -name *groups 2> /dev/null | parallel -j $(CPU) "grep -wf {} $$(find ${DIR}/orthofuse/${RUNOUT}/ -name contigs.csv 2> /dev/null) > {1}.orthout 2> /dev/null"
	echo About to delete all the text files
	find ${DIR}/orthofuse/${RUNOUT}/ -name *groups -delete
	echo Search output files
	find ${DIR}/orthofuse/${RUNOUT}/ -name *orthout 2> /dev/null | parallel -j $(CPU) "awk -F, 'BEGIN {max = 0} {if (\$$9>max) max=\$$9} END {print \$$1 \"\\t\" max}'" | tee -a ${DIR}/orthofuse/${RUNOUT}/good.list
	find ${DIR}/orthofuse/${RUNOUT}/ -name *orthout -delete
	python ${MAKEDIR}/scripts/filter.py ${DIR}/orthofuse/${RUNOUT}/merged.fasta <(awk '{print $$1}' ${DIR}/orthofuse/${RUNOUT}/good.list) > ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta
	cp ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	rm ${DIR}/orthofuse/${RUNOUT}/good.list ${DIR}/orthofuse/${RUNOUT}/merged.fasta

busco.done:
	python $$(which run_BUSCO.py) -i ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta -m transcriptome --cpu $(CPU) -l /mnt/lustre/macmaneslab/macmanes/BUSCODB/${LINEAGE} -o ${DIR}/reports/${RUNOUT} && \
	touch busco.done

transrate.done:
	cd ${DIR}/reports && \
	transrate -o transrate_${basename ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta .fasta}  -a ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta --left ${READ1} --right ${READ2} -t $(CPU) && \
	touch transrate.done

reportgen:
	printf "\n\n*****  QUALITY REPORT FOR: ${ASSEMBLY} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/run_${RUNOUT} -name short*) | sed -n 8p  | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	cat $$(find reports/transrate_${basename ${ASSEMBLY} .fasta} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${basename ${ASSEMBLY} .fasta}
	printf " \n\n"
