#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#       orthofuser.mk all READ1= READ2= CPU= RUNOUT= FASTADIR= LINEAGE=
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
BUSCO_CONFIG_FILE := ${MAKEDIR}/software/config.ini
export BUSCO_CONFIG_FILE



setup:${DIR}/ortho_setup.done
merge:${DIR}/orthofuse/${RUNOUT}/merged.fasta
orthotransrate:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
orthofusing:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
cdhit:${DIR}/assemblies/${RUNOUT}.ORP.fasta
busco:${DIR}/reports/${RUNOUT}.busco.done
transrate:${DIR}/reports/${RUNOUT}.transrate.done
salmon:${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf
clean:

all: setup merge orthotransrate orthofusing cdhit busco transrate salmon clean

.DELETE_ON_ERROR:
.PHONY:report

${DIR}/ortho_setup.done:
	@mkdir -p ${DIR}/reports
	@mkdir -p ${DIR}/quants
	touch ${DIR}/ortho_setup.done

${DIR}/orthofuse/${RUNOUT}/merged.fasta:
	mkdir -p ${DIR}/orthofuse/${RUNOUT}/working
	for fasta in $$(ls ${FASTADIR}); do python ${MAKEDIR}/scripts/long.seq.py ${FASTADIR}/$$fasta ${DIR}/orthofuse/${RUNOUT}/working/$$fasta.short.fasta 200; done
	( \
	source ${MAKEDIR}/software/anaconda/install/bin/activate py27; \
	python $$(which orthofuser.py) -I 4 -f ${DIR}/orthofuse/${RUNOUT}/working/ -og -t $(CPU) -a $(CPU); \
	source ${MAKEDIR}/software/anaconda/install/bin/activate orp_v2;\
	)
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
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' 2> /dev/null | parallel -j $(CPU) "awk -F, -v max=0 '{if(\$$14>max){want=\$$1; max=\$$14}}END{print want}'" >> ${DIR}/orthofuse/${RUNOUT}/good.list
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' -delete
	python ${MAKEDIR}/scripts/filter.py ${DIR}/orthofuse/${RUNOUT}/merged.fasta ${DIR}/orthofuse/${RUNOUT}/good.list > ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	rm ${DIR}/orthofuse/${RUNOUT}/good.list

${DIR}/assemblies/${RUNOUT}.ORP.fasta:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	cd ${DIR}/assemblies/ && cd-hit-est -M 5000 -T $(CPU) -c .98 -i ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta -o ${DIR}/assemblies/${RUNOUT}.ORP.fasta
	diamond blastx -p $(CPU) -e 1e-8 --top 0.1 -q ${DIR}/assemblies/${RUNOUT}.ORP.fasta -d ${MAKEDIR}/software/diamond/swissprot  -o ${DIR}/assemblies/${RUNOUT}.ORP.diamond.txt
	awk '{print $$2}' ${DIR}/assemblies/${RUNOUT}.ORP.diamond.txt | awk -F "|" '{print $$3}' | cut -d _ -f2 | sort | uniq | wc -l > ${DIR}/assemblies/${RUNOUT}.unique.ORP.txt
	rm ${DIR}/assemblies/${RUNOUT}.ORP.fasta.clstr

${DIR}/reports/${RUNOUT}.busco.done:${DIR}/assemblies/${RUNOUT}.ORP.fasta
	python $$(which run_BUSCO.py) -i ${FASTADIR}/${RUNOUT}.ORP.fasta -m transcriptome -f --cpu $(CPU) -o ${RUNOUT}
	mv run_${RUNOUT} ${DIR}/reports/
	touch ${DIR}/reports/${RUNOUT}busco.done

${DIR}/reports/${RUNOUT}.transrate.done:${DIR}/reports/${RUNOUT}.busco.done
	transrate -o ${DIR}/reports/transrate_${RUNOUT}  -a ${FASTADIR}/${RUNOUT}.ORP.fasta --left ${READ1} --right ${READ2} -t $(CPU)
	touch ${DIR}/reports/${RUNOUT}.transrate.done

${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf:${DIR}/reports/${RUNOUT}.transrate.done
	salmon index --no-version-check -t ${FASTADIR}/${RUNOUT}.ORP.fasta  -i ${RUNOUT}.ortho.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.ortho.idx --seqBias --gcBias -l a -1 ${READ1} -2 ${READ2} -o ${DIR}/quants/salmon_orthomerged_${RUNOUT}
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

	printf " \n\n"
	source ${MAKEDIR}/software/anaconda/install/bin/deactivate
