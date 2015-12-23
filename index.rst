==============================================
Oyster River Protocol For Transcriptome Assembly
==============================================

    The OR Protocol for transcriptome assembly is an actively developed, evidenced based method for optimizing transcriptome assembly. 

--------------------------------------------------
Contact Information for Professor Matthew MacManes
--------------------------------------------------

    - Gitter (preferred) https://gitter.im/macmanes-lab/Oyster_River_Protocol
    - Email (good): Matthew.MacManes@unh.edu
    - Twitter (good): @macmanes
    - Phone (discouraged): 603-862-4052
    - Office (I'm hiding under my desk): 189 Rudman Hall

--------------------------------------------------
 :doc:`aws_setup`
--------------------------------------------------

0. Archive Reads
-----------------------------------
::

  gzip *fastq

1. Initial Quality Check
-----------------------------------

::

  SolexaQA.pl ../right.fq
  

2. Error Correct
-----------------------------------

::

  perl run_rcorrector.pl -k 31 -t 30 \
  -1 ../reads/file_1.fastq \
  -2 ../reads/file_2.fastq

3. Assemble
-----------------------------------

::

  Trinity --seqType fq --max_memory 40G --trimmomatic --CPU 30\
  --left ../reads/file_1.cor.fastq --right ../reads/file_1.cor.fastq \
  --quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2   TRAILING:2 MINLEN:25"

4. Quality Check
-----------------------------------

5. Filter
-----------------------------------

6. Report
-----------------------------------
