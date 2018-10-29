SAMPLE=1
ASSEMBLER=shannon
for i in $(ls -lth reads/*1.fastq.gz); do
  F=$(basename $i _1.fastq.gz);
  $HOME/Oyster_River_Protocol/orp.mk report CPU=24 ASSEMBLY=$F.$SAMPLE.$ASSEMBLER.fasta \
  READ1=$F_1.fastq.gz READ2=$F_2.fastq.gz \
  LINEAGE=eukaryota SAMP=$SAMPLE --dry-run;
done




dataset=SRR3499127
for i in $(ls assemblies/$dataset.1.*fasta); do
  F=$(basename $i)
  $HOME/Oyster_River_Protocol/orp.mk report SAMP=1 CPU=24 ASSEMBLY=$F READ1=SRR3499127_1.fastq.gz READ2=SRR3499127_2.fastq.gz LINEAGE=eukaryota_odb9 --dry-run;
done
