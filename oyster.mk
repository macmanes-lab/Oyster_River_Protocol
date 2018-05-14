#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	oyster.mk main READ1= READ2= MEM=500 CPU=24 RUNOUT=runname
# oyster.mk orthofuse FASTADIR= READ1= READ2= MEM=500 CPU=24 RUNOUT=runname
#

VERSION = 1.2.0
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
shannonpath := $(shell which shannon.py 2>/dev/null)
brewpath := $(shell which brew 2>/dev/null)
rcorrpath := $(shell which rcorrector 2>/dev/null)
trimmomaticpath := $(shell which trimmomatic 2>/dev/null)
trinitypath := $(shell which Trinity 2>/dev/null)
spadespath := $(shell which rnaspades.py 2>/dev/null)
salmonpath := $(shell which salmon 2>/dev/null)
mclpath := $(shell which mcl 2>/dev/null)
buscopath := $(shell which run_BUSCO.py 2>/dev/null)
seqtkpath := $(shell which seqtk 2>/dev/null)
transratepath := $(shell which transrate 2>/dev/null)




run_trimmomatic:
run_rcorrector:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
run_trinity:${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta
run_spades55:${DIR}/assemblies/${RUNOUT}.spades55.fasta
run_spades75:${DIR}/assemblies/${RUNOUT}.spades75.fasta
run_shannon:${DIR}/assemblies/${RUNOUT}.shannon.fasta
merge:${DIR}/orthofuse/${RUNOUT}/merged.fasta
orthotransrate:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
orthofusing:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
salmon:${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf
shmlast:${DIR}/assemblies/shmlast/${RUNOUT}.trinity.crbl.csv
posthack:${DIR}/assemblies/shmlast/newbies.fasta
main: setup check welcome run_trimmomatic run_rcorrector run_trinity run_spades75 run_spades55 run_shannon merge orthotransrate orthofusing shmlast posthack salmon busco transrate report
orthofuse:merge orthotransrate orthofusing
report:busco transrate reportgen
busco:${DIR}/reports/busco.done
transrate:${DIR}/reports/transrate.done

.DELETE_ON_ERROR:
.PHONY:report check

setup:
	@mkdir -p ${DIR}/scripts
	@mkdir -p ${DIR}/reads
	@mkdir -p ${DIR}/assemblies
	@mkdir -p ${DIR}/rcorr
	@mkdir -p ${DIR}/reports
	@mkdir -p ${DIR}/orthofuse
	@mkdir -p ${DIR}/quants

check:
ifdef salmonpath
else
	$error("*** SALMON is not installed, must fix ***")
endif
ifdef transratepath
else
	$error("*** TRANSRATE is not installed, must fix ***")
endif
ifdef seqtkpath
else
	$error("*** SEQTK is not installed, must fix ***")
endif
ifdef buscopath
else
	$error("*** BUSCO is not installed, must fix ***")
endif
ifdef mclpath
else
	$error("*** MCL is not installed, must fix ***")
endif
ifdef spadespath
else
	$error("*** SPADES is not installed, must fix ***")
endif
ifdef trinitypath
else
	$error("*** TRINITY is not installed, must fix ***")
endif
ifdef trimmomaticpath
else
	@echo "Maybe TRIMMOMATIC is not installed, or maybe you are working on Bridges"
endif
ifdef shannonpath
else
	$error("*** SHANNON is not installed, must fix ***")
endif
ifdef rcorrpath
else
	$error("*** RCORRECTOR is not installed, must fix ***")
endif


welcome:
	printf "\n\n*****  Welcome to the Oyster River ***** \n"
	printf "*****  This is version ${VERSION} ***** \n\n "
	printf " \n\n"

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
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta -i ${RUNOUT}.salmon.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.salmon.idx --seqBias --gcBias -l a -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -o ${DIR}/quants/salmon_trin_${RUNOUT}
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta -m transcriptome --cpu $(CPU) -o ${RUNOUT}.trin
	rm -fr ${RUNOUT}.salmon.idx

${DIR}/assemblies/${RUNOUT}.spades55.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	rnaspades.py --only-assembler -o ${DIR}/assemblies/${RUNOUT}.spades_k55 --threads $(CPU) --memory $(MEM) -k 55 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq
	mv ${DIR}/assemblies/${RUNOUT}.spades_k55/transcripts.fasta ${DIR}/assemblies/${RUNOUT}.spades55.fasta
	rm -fr ${DIR}/assemblies/${RUNOUT}.spades_k55
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.spades55.fasta -i ${RUNOUT}.salmon.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.salmon.idx --seqBias --gcBias -l a -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -o ${DIR}/quants/salmon_sp55_${RUNOUT}
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${RUNOUT}.spades55.fasta -m transcriptome --cpu $(CPU) -o ${RUNOUT}.spades55
	rm -fr ${RUNOUT}.salmon.idx


${DIR}/assemblies/${RUNOUT}.spades75.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	rnaspades.py --only-assembler -o ${DIR}/assemblies/${RUNOUT}.spades_k75 --threads $(CPU) --memory $(MEM) -k 75 -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq
	mv ${DIR}/assemblies/${RUNOUT}.spades_k75/transcripts.fasta ${DIR}/assemblies/${RUNOUT}.spades75.fasta
	rm -fr ${DIR}/assemblies/${RUNOUT}.spades_k75
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.spades75.fasta -i ${RUNOUT}.salmon.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.salmon.idx --seqBias --gcBias -l a -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -o ${DIR}/quants/salmon_sp75_${RUNOUT}
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${RUNOUT}.spades75.fasta -m transcriptome --cpu $(CPU) -o ${RUNOUT}.spades75
	rm -fr ${RUNOUT}.salmon.idx


${DIR}/assemblies/${RUNOUT}.shannon.fasta:${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq
	rm -fr ${DIR}/assemblies/${RUNOUT}.shannon
	python $$(which shannon.py) -o ${DIR}/assemblies/${RUNOUT}.shannon --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -p $(CPU) -K 75
	mv ${DIR}/assemblies/${RUNOUT}.shannon/shannon.fasta ${DIR}/assemblies/${RUNOUT}.shannon.fasta
	rm -fr ${DIR}/assemblies/${RUNOUT}.shannon
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.shannon.fasta -i ${RUNOUT}.salmon.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.salmon.idx --seqBias --gcBias -l a -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -o ${DIR}/quants/salmon_shannon_${RUNOUT}
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${RUNOUT}.shannon.fasta -m transcriptome --cpu $(CPU) -o ${RUNOUT}.shannon
	rm -fr ${RUNOUT}.salmon.idx


${DIR}/orthofuse/${RUNOUT}/merged.fasta:
	mkdir -p ${DIR}/orthofuse/${RUNOUT}/working
	for fasta in $$(ls ${DIR}/assemblies/${RUNOUT}*fasta); do python ${MAKEDIR}/scripts/long.seq.py ${DIR}/assemblies/$$(basename $$fasta) ${DIR}/orthofuse/${RUNOUT}/working/$$(basename $$fasta).short.fasta 200; done
	python $$(which orthofuser.py) -I 4 -f ${DIR}/orthofuse/${RUNOUT}/working/ -og -t $(CPU) -a $(CPU)
	cat ${DIR}/orthofuse/${RUNOUT}/working/*short.fasta > ${DIR}/orthofuse/${RUNOUT}/merged.fasta

${DIR}/orthofuse/${RUNOUT}/orthotransrate.done:${DIR}/orthofuse/${RUNOUT}/merged.fasta
	export END=$$(wc -l $$(find ${DIR}/orthofuse/${RUNOUT}/working/ -name Orthogroups.txt 2> /dev/null) | awk '{print $$1}') && \
	export ORTHOINPUT=$$(find ${DIR}/orthofuse/${RUNOUT}/working/ -name Orthogroups.txt 2> /dev/null) && \
	echo $$(eval echo "{1..$$END}") | tr ' ' '\n' > ${DIR}/orthofuse/${RUNOUT}/list && \
	cat ${DIR}/orthofuse/${RUNOUT}/list | parallel  -j $(CPU) -k "sed -n ''{}'p' $$ORTHOINPUT | tr ' ' '\n' | sed '1d' > ${DIR}/orthofuse/${RUNOUT}/{1}.groups"
	${MAKEDIR}/software/orp-transrate/transrate -o ${DIR}/orthofuse/${RUNOUT}/merged -t $(CPU) -a ${DIR}/orthofuse/${RUNOUT}/merged.fasta --left ${READ1} --right ${READ2}
	touch ${DIR}/orthofuse/${RUNOUT}/orthotransrate.done

${DIR}/assemblies/${RUNOUT}.orthomerged.fasta:${DIR}/orthofuse/${RUNOUT}/orthotransrate.done
	echo All the text files are made, start GREP
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*groups' 2> /dev/null | parallel -j $(CPU) "grep -wf {} $$(find ${DIR}/orthofuse/${RUNOUT}/ -name contigs.csv 2> /dev/null) > {1}.orthout 2> /dev/null"
	echo About to delete all the text files
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*groups' -delete
	echo Search output files
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' 2> /dev/null | parallel -j $(CPU) "awk -F, -v max=0 '{if(\$$14>max){want=\$$1; max=\$$14}}END{print want}'" >> ${DIR}/orthofuse/${RUNOUT}/good.list
	find ${DIR}/orthofuse/${RUNOUT}/ -name '*orthout' -delete
	python ${MAKEDIR}/scripts/filter.py ${DIR}/orthofuse/${RUNOUT}/merged.fasta ${DIR}/orthofuse/${RUNOUT}/good.list > ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta
	cp ${DIR}/orthofuse/${RUNOUT}/${RUNOUT}.orthomerged.fasta ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	rm ${DIR}/orthofuse/${RUNOUT}/good.list

${DIR}/assemblies/shmlast/${RUNOUT}.trinity.crbl.csv:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta ${DIR}/assemblies/${RUNOUT}.shannon.fasta ${DIR}/assemblies/${RUNOUT}.spades75.fasta ${DIR}/assemblies/${RUNOUT}.spades55.fasta ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta
	mkdir -p ${DIR}/assemblies/shmlast
	( \
	source ${MAKEDIR}/software/anaconda/install/bin/activate; \
	shmlast crbl -q ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta -d ${MAKEDIR}/software/shmlast/uniprot_sprot.fasta --n_threads $(CPU) -o ${DIR}/assemblies/shmlast/${RUNOUT}.orthomerged.crbl.csv ; \
	shmlast crbl -q ${DIR}/assemblies/${RUNOUT}.shannon.fasta -d ${MAKEDIR}/software/shmlast/uniprot_sprot.fasta --n_threads $(CPU) -o ${DIR}/assemblies/shmlast/${RUNOUT}.shannon.crbl.csv; \
	shmlast crbl -q ${DIR}/assemblies/${RUNOUT}.spades75.fasta -d ${MAKEDIR}/software/shmlast/uniprot_sprot.fasta --n_threads $(CPU) -o ${DIR}/assemblies/shmlast/${RUNOUT}.spades75.crbl.csv; \
	shmlast crbl -q ${DIR}/assemblies/${RUNOUT}.spades55.fasta -d ${MAKEDIR}/software/shmlast/uniprot_sprot.fasta --n_threads $(CPU) -o ${DIR}/assemblies/shmlast/${RUNOUT}.spades55.crbl.csv; \
	shmlast crbl -q ${DIR}/assemblies/${RUNOUT}.trinity.Trinity.fasta -d ${MAKEDIR}/software/shmlast/uniprot_sprot.fasta --n_threads $(CPU) -o ${DIR}/assemblies/shmlast/${RUNOUT}.trinity.crbl.csv; \
	source deactivate; \
	)

${DIR}/assemblies/shmlast/newbies.fasta:${DIR}/assemblies/shmlast/${RUNOUT}.trinity.crbl.csv
	cd ${DIR}/assemblies/shmlast/ && cut -d, -f14 ${RUNOUT}.orthomerged.crbl.csv | cut -d "|" -f3 | cut -d "_" -f1 | sort --parallel=20 |uniq > list1
	cd ${DIR}/assemblies/shmlast/ && cut -d, -f14 ${RUNOUT}.{shannon,spades75,spades55,trinity}.crbl.csv | cut -d "|" -f3 | cut -d "_" -f1 | sort --parallel=20 |uniq > list2
	cd ${DIR}/assemblies/shmlast/ && grep -vwf list1 list2 > list3
	cd ${DIR}/assemblies/shmlast/ && for item in $$(cat list3); do grep $item ${RUNOUT}.{shannon,spades75,spades55,trinity}.crbl.csv | head -1 | cut -d, -f9 >> list4 ; done; sort list4 | uniq > list5
	cd ${DIR}/assemblies/shmlast/ && python ${MAKEDIR}/scripts/filter.py <(cat ../${RUNOUT}.{spades55,spades75,shannon,trinity.Trinity}.fasta) list5 >> newbies.fasta


${DIR}/reports/busco.done:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	python $$(which run_BUSCO.py) -i ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta -m transcriptome --cpu $(CPU) -o ${RUNOUT}.orthomerged
	mv run_${RUNOUT}* ${DIR}/reports/
	touch ${DIR}/reports/${RUNOUT}.busco.done

${DIR}/reports/transrate.done:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	${MAKEDIR}/software/orp-transrate/transrate -o ${DIR}/reports/transrate_${RUNOUT}  -a ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta --left ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq --right ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -t $(CPU)
	touch ${DIR}/reports/${RUNOUT}.transrate.done

${DIR}/quants/salmon_orthomerged_${RUNOUT}/quant.sf:${DIR}/assemblies/${RUNOUT}.orthomerged.fasta
	salmon index --no-version-check -t ${DIR}/assemblies/${RUNOUT}.orthomerged.fasta  -i ${RUNOUT}.ortho.idx --type quasi -k 31
	salmon quant --no-version-check -p $(CPU) -i ${RUNOUT}.ortho.idx --seqBias --gcBias -l a -1 ${DIR}/rcorr/${RUNOUT}.TRIM_1P.cor.fq -2 ${DIR}/rcorr/${RUNOUT}.TRIM_2P.cor.fq -o ${DIR}/quants/salmon_orthomerged_${RUNOUT}
	rm -fr ${RUNOUT}.ortho.idx

reportgen:
	printf "\n\n*****  QUALITY REPORT FOR: ${RUNOUT} using the ORP version ${VERSION} **** \n\n"
	printf "*****  BUSCO SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/run_${RUNOUT}.orthomerged -name 'short*') | sed -n 8p  | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE SCORE ~~~~~>           " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$37}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf "*****  TRANSRATE OPTIMAL SCORE ~~~~~>   " | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	cat $$(find reports/transrate_${RUNOUT} -name assemblies.csv) | awk -F , '{print $$38}' | sed -n 2p | tee -a ${DIR}/reports/qualreport.${RUNOUT}
	printf " \n\n"
