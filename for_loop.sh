SAMPLE=1
ASSEMBLER=shannon
for i in $(ls -lth reads/*1.fastq.gz); do
  F=$(basename $i _1.fastq.gz);
  $HOME/Oyster_River_Protocol/orp.mk report CPU=24 ASSEMBLY=$F.$SAMPLE.$ASSEMBLER.fasta \
  READ1=$F_1.fastq.gz READ2=$F_2.fastq.gz \
  LINEAGE=eukaryota SAMP=$SAMPLE --dry-run;
done
