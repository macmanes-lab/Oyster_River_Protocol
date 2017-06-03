input=Orthogroups.txt
END=$(wc -l $input | awk '{print $1}')
START=1

/share/orthonucl/OrthoFinder/orthofinder/orthofinder.py -f fasta/ -og -t 40 -a 40



for i in $(eval echo "{$START..$END}") ; do
  sed -n ''$i'p' $input | tr ' ' '\n' | grep -f - trans/merged/contigs.csv \
  | awk -F, '{print $1 "\t" $9}' | sort -rnk2,2 | head -1 | tee -a good.list
done
