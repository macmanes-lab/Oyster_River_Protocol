#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#       orthofuser.mk all READ1= READ2= CPU= RUNOUT= FASTADIR= LINEAGE= TPM_FILT=
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
VERSION := ${shell cat  ${MAKEDIR}version.txt}
TPM_FILT =
INPUT_FASTAS := ${shell ls ${DIR}/${FASTADIR}}


setup:${DIR}/ortho_setup.done
merge:${DIR}/orthofuse/${RUNOUT}/merged.fasta
orthotransrate:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
orthofusing:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
diamond:${DIR}/assemblies/diamond/diamond.done
posthack:${DIR}/assemblies/diamond/${RUNOUT}.newbies.fasta
cdhit:${DIR}/assemblies/${RUNOUT}.ORP.fasta
filter:${DIR}/assemblies/${RUNOUT}.filter.done
busco:${DIR}/reports/${RUNOUT}.busco.done
transrate:${DIR}/reports/${RUNOUT}.transrate.done
salmon:${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf
reportgen:


all: setup merge orthotransrate orthofusing diamond posthack cdhit salmon filter busco transrate reportgen

.DELETE_ON_ERROR:
.PHONY:report

${DIR}/ortho_setup.done:
	@mkdir -p ${DIR}/reports
	@mkdir -p ${DIR}/quants
	@mkdir -p ${DIR}/assemblies
	touch ${DIR}/ortho_setup.done
	@mkdir -p ${DIR}/assemblies/working
	@mkdir -p ${DIR}/assemblies/diamond

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

${DIR}/assemblies/diamond/diamond.done:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
        echo Starting diamond
        for fasta in $$(ls ${FASTADIR}); do echo Running diamond on $$fasta assembly; diamond blastx -p $(CPU) -e 1e-8 --top 0.1 -q ${DIR}/${FASTADIR}/$$fasta -d ${MAKEDIR}/software/diamond/swissprot -o ${DIR}/assemblies/diamond/$$(basename $$fasta).inputfasta.diamond.txt; awk '{print $$2}' ${DIR}/assemblies/diamond/$$(basename $$fasta).inputfasta.diamond.txt | awk -F "|" '{print $$3}' | cut -d _ -f2 | sort | uniq | wc -l > ${DIR}/assemblies/diamond/$$(basename $$fasta).unique.txt; done
	# run for orthomerged assembly 
	diamond blastx -p $(CPU) -e 1e-8 --top 0.1 -q ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta -d ${MAKEDIR}/software/diamond/swissprot -o ${DIR}/assemblies/diamond/${RUNOUT}.orthomerged.diamond.txt)
	awk '{print $$2}' ${DIR}/assemblies/diamond/${RUNOUT}.orthomerged.diamond.txt | awk -F "|" '{print $$3}' | cut -d _ -f2 | sort | uniq | wc -l > ${DIR}/assemblies/diamond/${RUNOUT}.orthomerged.unique.txt
	# touch to signal end of command.
	touch ${DIR}/assemblies/diamond/diamond.done

# list1 is unique geneIDs in orthomerged
#list2 is unique  geneIDs in other assemblies
#list3 is which genes are in other assemblies but not in orthomerged
#list4 are the contig IDs from input_fastas from list3
#list5 is unique IDs contig IDs from list4
#list6 is contig IDs from orthomerged FASTAs
#list7 is stuff that is in 5 but not 6

${DIR}/assemblies/diamond/${RUNOUT}.newbies.fasta:${DIR}/assemblies/diamond/diamond.done
	cd ${DIR}/assemblies/diamond/ && cut -f2 ${RUNOUT}.orthomerged.diamond.txt | cut -d "|" -f3 | cut -d "_" -f1 | sort --parallel=20 |uniq > ${RUNOUT}.list1
	cd ${DIR}/assemblies/diamond/ && touch inputfastaIDs.diamond.txt
	cd ${DIR}/assemblies/diamond/ && for file in $$(ls *inputfasta.diamond.txt); do cat $$file >> inputfastaIDs.diamond.txt; done
	cd ${DIR}/assemblies/diamond/ && cut -f2 inputfastaIDs.diamond.txt | cut -d "|" -f3 | cut -d "_" -f1 | sort --parallel=20 | uniq > ${RUNOUT}.list2
	cd ${DIR}/assemblies/diamond/ && grep -xFvwf ${RUNOUT}.list1 ${RUNOUT}.list2 > ${RUNOUT}.list3
	cd ${DIR}/assemblies/diamond/ && for item in $$(cat ${RUNOUT}.list3); do grep -F $$item inputfastaIDs.diamond.txt | head -1 | cut -f1; done | cut -d ":" -f2 | sort | uniq >> ${RUNOUT}.list5
	cd ${DIR}/assemblies/diamond/ && grep -F ">" ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta | sed 's_>__' > ${RUNOUT}.list6
	cd ${DIR}/assemblies/diamond/ && grep -xFvwf ${RUNOUT}.list6 ${RUNOUT}.list5 > ${RUNOUT}.list7
	cd ${DIR}/assemblies/diamond/ && python ${MAKEDIR}/scripts/filter.py <(cat ${DIR}/${FASTADIR}/*) ${RUNOUT}.list7 >> ${RUNOUT}.newbies.fasta
	cd ${DIR}/assemblies/diamond/ &&  cat ${RUNOUT}.newbies.fasta ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta > tmp.fasta && mv tmp.fasta ${DIR}/assemblies/working/${RUNOUT}.orthomerged.fasta
	cd ${DIR}/assemblies/diamond/ && rm -f ${RUNOUT}.list*

${DIR}/assemblies/${RUNOUT}.ORP.fasta:${DIR}/assemblies/diamond/${RUNOUT}.newbies.fasta
	cd ${DIR}/assemblies/ && cd-hit-est -M 5000 -T $(CPU) -c .98 -i ${DIR}/assemblies/working/${RUNOUT}.orthomerged.fasta -o ${DIR}/assemblies/${RUNOUT}.ORP.fasta
	diamond blastx -p $(CPU) -e 1e-8 --top 0.1 -q ${DIR}/assemblies/${RUNOUT}.ORP.fasta -d ${MAKEDIR}/software/diamond/swissprot  -o ${DIR}/assemblies/${RUNOUT}.ORP.diamond.txt
	awk '{print $$2}' ${DIR}/assemblies/${RUNOUT}.ORP.diamond.txt | awk -F "|" '{print $$3}' | cut -d _ -f2 | sort | uniq | wc -l > ${DIR}/assemblies/working/${RUNOUT}.unique.ORP.txt
	rm ${DIR}/assemblies/${RUNOUT}.ORP.fasta.clstr
	
${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf:${DIR}/assemblies/${RUNOUT}.ORP.fasta
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.ORP.fasta  -i ${RUNOUT}.ortho.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.ortho.idx --seqBias --gcBias -l a -1 ${READ1} -2 ${READ2} -o ${DIR}/quants/salmon_orthomerged_${RUNOUT}
	rm -fr ${RUNOUT}.ortho.idx

${DIR}/assemblies/${RUNOUT}.filter.done ${DIR}/assemblies/working/${RUNOUT}.saveme.fasta:${DIR}/assemblies/${RUNOUT}.ORP.fasta ${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf
ifdef TPM_FILT
	echo Filtering out transcripts with lower than $TPM_FILT TPM
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

${DIR}/reports/${RUNOUT}.busco.done:${DIR}/assemblies/${RUNOUT}.filter.done
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${RUNOUT}.ORP.fasta -m transcriptome -f --cpu $(CPU) -o ${RUNOUT}
	mv run_${RUNOUT} ${DIR}/reports/
	touch ${DIR}/reports/${RUNOUT}.busco.done

${DIR}/reports/${RUNOUT}.transrate.done:${DIR}/reports/${RUNOUT}.busco.done
	transrate -o ${DIR}/reports/transrate_${RUNOUT}  -a ${DIR}/assemblies/${RUNOUT}.ORP.fasta --left ${READ1} --right ${READ2} -t $(CPU)
	touch ${DIR}/reports/${RUNOUT}.transrate.done


reportgen:
	printf "\n\n*****  QUALITY REPORT FOR: ${RUNOUT} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/run_${RUNOUT} -name 'short*') | sed -n 8p  | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf " \n\n"

	printf " \n Orthofuser complete \n"
	source ${MAKEDIR}/software/anaconda/install/bin/deactivate
