### Kmers as an assembly metric

Basically, can we learn something about reads, by quantitating the number of matches to a refseq database?

Here are the putative steps.

1. make fasta file of kmers from refseq. The specific data base is based on the taxa unser study. I'll use the mammalian RefSeq for starters.

2. Make jellyfish database of kmers (same length as RefSeq database)

3. Query the overlap. As the read set improves, the degree of overlap increases. This degree should plateau at some specific sequencing depth.

4. At the "right" kmer length, there should be close to 0 kmers found in the invert RefSeq database 

```
#!/bin/bash

cd /mnt/lustre/macmaneslab/macmanes/kmers



SUBSAMP=10000
SAMP=10k
READS=reads.fq
TAXA=vert
DATABASE="$TAXA".rna.fna

seqtk sample -s23 SRR2086412.TRIM_1P.cor.fq "$SUBSAMP" > "$SAMP"_"$TAXA"_"$READS"

for kmer in `seq 31 10 41`; do
        jellyfish count -m "$kmer" -s 100M -t 24 -o /dev/stdout -C "$DATABASE" | jellyfish dump /dev/stdin > "$kmer"mers_in_"$TAXA"_refseq.fasta;
        jellyfish count -m "$kmer" -s 100M -t 24 -o k"$kmer"_"$SAMP"reads.jf -C "$SAMP"_"$TAXA"_"$READS";
        jellyfish query k"$kmer"_"$SAMP"reads.jf -s "$kmer"mers_in_"$TAXA"_refseq.fasta | awk '$2 > 0' | tee -a k"$kmer"_"$SAMP"_matches.list;
done
```
