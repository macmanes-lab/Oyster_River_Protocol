# General usage for orthofuser.mk

This make file is meant to merge multiple transcriptomes together. In practice, this is the same functionality as when assemblies from Trinity, SPades, and TransAbyss are merged together in  `oyster.mk`. However, this make file was created specifically to increase the general applicability of assembled transcriptomes, and allow a user to merge together transcriptomes that have already been assembled (for example, from multiple species or treatments).

# Usage:
First, the user needs to activate the relevant conda environment using `source activate orp_v2`

Then the reads are specified with `READ1=` and `READ2=`. We suggest that the user uses files that contain the readsets used to produce all of the transcriptomes merged into a single file. The number of threads are specified using `CPU=`, obviously you do not want to specify more threads than you have access to (and if you are using a personal machine that you plan to continue working on you will want to specify fewer than the machine has or you are going to be frustrated). Similarly, you will want to b careful of the amount of memory you allow via `MEM=`. Giving too much will slow down the program, and might get an automatic kill if you are using a manager like SLURM. The `FASTADIR=` argument specifies the directory of the transcriptomes you want to merge reside. Importantly, `orthofuser.mk` will merge ALL the fasta files within that folder so be careful about that. The `RUNOUT=` argument specifies the naming scheme/output.

An example command to be executed:

```
source activate orp_v2

/PATH_TO_ORP/Oyster_River_Protocol/orthofuser.mk all \
READ1=merged_reads_1P.fq READ2=merged_reads_2P.fq CPU=24 MEM=500 \
RUNOUT=multispecies FASTADIR=assemblies \
LINEAGE=/PATH_TO_ORP/Oyster_River_Protocol/busco_dbs/eukaryota_odb9
```

# An important consideration:
This script will fail if you have multiple sequences with the same header. This is common when you are producing multiple assemblies via the same pipeline/assembler. You will need to change those prior to invoking `orthofuser.mk` or it will fail.

Some basic code to change the header of a fasta, assuming that your fasta is `YourAmazingTranscriptomeFromORP.fasta`:
```
awk '/^>/{print ">Transcript_" ++i; next}{print}' YourAmazingTranscriptomeFromORP.fasta > tmp.fasta
mv tmp.fasta > YourAmazingTranscriptomeFromORP.fasta
rm tmp.fasta
```
