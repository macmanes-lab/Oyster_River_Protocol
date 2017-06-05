input=Orthogroups.txt
END=$(wc -l $input | awk '{print $1}')
START=1

/share/orthonucl/OrthoFinder/orthofinder/orthofinder.py -f fasta/ -og -t 40 -a 40

cat fasta/*fasta > merged.fasta

/share/test/share/transrate-1.0.4beta/transrate -o merged -t 48 \
-a merged.fasta \
--left ../../SRR1199463_1.fastq.gz  \
--right ../../SRR1199463_2.fastq.gz

for i in $(eval echo "{$START..$END}") ; do
  sed -n ''$i'p' $input | tr ' ' '\n' | grep -f - trans/merged/contigs.csv \
  | awk -F, '{print $1 "\t" $9}' | sort -rnk2,2 | head -1 | tee -a good.list
done

python /mnt/data3/macmanes/Mc_Transcriptome/final/Sept16/filter.py merged.fasta <(awk '{print $1}' good.list)  >  orthomerged.fasta

python3 /share/busco/busco/scripts/run_BUSCO.py -i orthomerged.fasta -o mer -m tran -l $RAID/busco_dbs/eukaryota_odb9 -c 48

/share/test/share/transrate-1.0.4beta/transrate -o orthomerged -t 48 \
-a orthomerged.fasta \
--left ../../SRR1199463_1.fastq.gz  \
--right ../../SRR1199463_2.fastq.gz
