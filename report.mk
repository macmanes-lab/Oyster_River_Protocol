#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#  $HOME/Oyster_River/Protocol/report.mk main CPU=24 \
#  ASSEMBLY=test.fasta \
#  READ1=1.subsamp_1.cor.fq \
#  READ2=1.subsamp_2.cor.fq \
#  LINEAGE=eukaryota_odb9 \
#  RUNOUT=test
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=24
MEM=120
READ1=
READ2=
BUSCO := ${shell which run_BUSCO.py}
BUSCODIR := $(dir $(firstword $(BUSCO)))
RUNOUT =
ASSEMBLY=
LINEAGE=
BUSCODBDIR := ${MAKEDIR}/busco_dbs/
BUSCOUT := BUSCO_$(shell basename ${ASSEMBLY} .fasta)
salmonpath := $(shell which salmon 2>/dev/null)
buscopath := $(shell which run_BUSCO.py 2>/dev/null)
transratepath := $(shell which transrate 2>/dev/null)
BUSCO_CONFIG_FILE := ${MAKEDIR}/software/config.ini
export BUSCO_CONFIG_FILE
VERSION := ${shell cat  ${MAKEDIR}version.txt}

help:
main: setup check welcome diamond busco transrate strandeval reportgen
diamond:${DIR}/reports/${RUNOUT}.unique.txt
busco:${DIR}/reports/${RUNOUT}.busco.done
transrate:${DIR}/reports/${RUNOUT}.transrate.done
clean:
setup:${DIR}/setup.done
strandeval:{DIR}/reports/${RUNOUT}.strandeval.done

.DELETE_ON_ERROR:
.PHONY:report check clean

${DIR}/setup.done:
	@mkdir -p ${DIR}/reports
	touch ${DIR}/setup.done

check:
ifdef salmonpath
else
	$(error "\n\n*** SALMON is not installed, must fix ***")
endif
ifdef transratepath
else
	$(error "\n\n*** TRANSRATE is not installed, must fix ***")
endif
ifdef buscopath
else
	$(error "\n\n*** BUSCO is not installed, must fix ***")
endif


help:
	printf "\n\n*****  Welcome to the Oyster River Report Generation Tool ***** \n"
	printf "*****  This is version ${VERSION} *****\n\n"
	printf "Usage:\n\n"
	printf "/path/to/Oyster_River/Protocol/report.mk main CPU=24\n"
	printf "ASSEMBLY=test.fasta\n"
	printf "LINEAGE=eukaryota_odb9\n"
	printf "READ1=1.subsamp_1.cor.fq\n"
	printf "READ2=1.subsamp_2.cor.fq\n"
	printf "RUNOUT=test\n\n"

welcome:
	printf "\n\n*****  Welcome to the Oyster River Report Generation Tool ***** \n"
	printf "*****  This is version ${VERSION} ***** \n\n "
	printf " \n\n"


${DIR}/reports/${RUNOUT}.busco.done:${ASSEMBLY}
	export BUSCO_CONFIG_FILE=${MAKEDIR}/software/config.ini
	python $$(which run_BUSCO.py) -i ${ASSEMBLY} -m transcriptome --cpu $(CPU) -o ${RUNOUT} --lineage_path ${BUSCODBDIR}/${LINEAGE}
	mv run_${RUNOUT} ${DIR}/reports/
	touch ${DIR}/reports/${RUNOUT}.busco.done

${DIR}/reports/${RUNOUT}.transrate.done:${ASSEMBLY}
	${MAKEDIR}/software/orp-transrate/transrate -o ${DIR}/reports/transrate_${RUNOUT}  -a ${ASSEMBLY} --left ${READ1} --right ${READ2} -t $(CPU)
	touch ${DIR}/reports/${RUNOUT}.transrate.done
	find ${DIR}/reports/transrate_${RUNOUT}/ -name "*bam" -delete

${DIR}/reports/${RUNOUT}.unique.txt:${ASSEMBLY}
	diamond blastx -p $(CPU) -e 1e-8 --top 0.1 -q ${ASSEMBLY} -d ${MAKEDIR}/software/diamond/swissprot -o ${DIR}/reports/${RUNOUT}.diamond.txt
	awk '{print $$2}' ${DIR}/reports/${RUNOUT}.diamond.txt | awk -F "|" '{print $$3}' | cut -d _ -f2 | sort | uniq | wc -l > ${DIR}/reports/${RUNOUT}.unique.txt


clean:
	rm -fr ${DIR}/reports/busco.done ${DIR}/reports/transrate.done ${DIR}/reports/${RUNOUT}.unique.txt ${DIR}/reports/run_${RUNOUT} ${DIR}/reports/transrate_${RUNOUT}/

{DIR}/reports/${RUNOUT}.strandeval.done:
	bwa index -p ${RUNOUT} ${ASSEMBLY}
	bwa mem -t $(CPU) ${RUNOUT} \
	<(seqtk sample -s 23894 ${READ1} 200000) \
	<(seqtk sample -s 23894 ${READ2} 200000) \
	| samtools view -@10 -Sb - \
	| samtools sort -T ${RUNOUT} -O bam -@10 -o "${RUNOUT}".sorted.bam -
	perl -I $$(dirname $$(readlink -f $$(which Trinity)))/PerlLib ${MAKEDIR}/scripts/examine_strand.pl "${RUNOUT}".sorted.bam ${RUNOUT}
	hist  -p '#' -c red <(cat ${RUNOUT}.dat | awk '{print $$5}' | sed  1d)
	rm -f "${RUNOUT}".sorted.bam
	touch ${DIR}/reports/${RUNOUT}.strandeval.done
	printf "\n\n*****  See the following link for interpretation ***** \n"
	printf "*****  LINK ***** \n\n"

reportgen:
	printf "\n\n*****  QUALITY REPORT FOR: ${RUNOUT} ****"
	printf "\n*****  THE ASSEMBLY CAN BE FOUND HERE: ${ASSEMBLY} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find ${DIR}/reports/run_${RUNOUT} -name 'short*') | sed -n 8p  | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find ${DIR}/reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find ${DIR}/reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  UNIQUE GENES ~~~~~~~~~>          " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat ${DIR}/reports/${RUNOUT}.unique.txt | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf " \n\n"
