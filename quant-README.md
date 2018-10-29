# General usage for quant.mk

This make file is meant to be used to trim adaptors, error correct reads, and pseudo-map to a reference transcriptome (for example one produced by the Oyster River Protocol...) for an individul sample. The user specifies some sample information, then the make file runs through to quantification of that sample. 

# Usage:
First, the user needs to activate the relevant conda environment using `source activate orp_v2`

Then the sample is specified with `SAMPLE=` and the ending of the file name with `SUFFIX=`. For example, a paired-end sample might contain two files with the forward reads named `Brie_1P.fq.gz` and the reverse reads named `Brie_2P.fq.gz`. In this case, the user would specify `SAMPLE=Brie_` and `SUFFIX=P.fq.gz`. Basically the sample is the unique identifier up to the indication of the forward or reverse read, and the suffix is everything after that designation. This will work with either compressed or uncompressed fastqs. The reference transcriptome is specified with `TRANSCRIPTOME=`

The number of threads are specified using `CPU=`, obviously you do not want to specify more threads than you have access to (and if you are using a personal machine that you plan to continue working on you will want to specify fewer than the machine has or you are going to be frustrated). Similarly, you will want to b careful of the amount of memory you allow via `MEM=`. Giving too much will slow down the program, and might get an automatic kill if you are using a manager like SLURM. 

An example command to be executed:

```
source activate orp_v2

/PATH_TO_ORP/Oyster_River_Protocol/quant.mk all \
MEM=500 CPU=24 \
SAMPLE=Brie_ \
SUFFIX=.fq.gz \
TRANSCRIPTOME=DeliciousCheeseFungi.fa
```

# Running ALL THE SAMPLES:
This make file allows the user to easily run through many samples/replicates from an experiment concurrently. You can easily supply a number of read files to map to a transcriptome with a simple for loop. For example, you have the following samples to map to a transcriptome:

```
Brie_1P.fq.gz
Brie_2P.fq.gz
Gouda_1P.fq.gz
Gouda_2P.fq.gz
Parmesan_1P.fq.gz
Parmesan_2P.fq.gz
```

You could set the individual samples into a variable (manually or programmatically) and feed them in to `quant.mk` in a for loop:

```
cheeses=$(Brie Gouda Parmesan)

source activate orp_v2

for i in $cheeses
do
/PATH_TO_ORP/Oyster_River_Protocol/quant.mk all \
MEM=500 CPU=24 \
SAMPLE=$i \
SUFFIX=.fq.gz \
TRANSCRIPTOME=DeliciousCheeseFungi.fa
done
```

# An important consideration:
This script will give you erroneous output if you do not correctly specify your adaptors in the `barcodes.fa` file in the `PATH_TO_ORP/Oyster_River_Protocol/barcodes/` directory. Presumably this will not be an issue here because the user will have already done this prior to assembling the initial transcriptome using the Oyster River Protocol.
