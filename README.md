# Oyster River Protocol

[![Join the chat at https://gitter.im/macmanes-lab/Oyster_River_Protocol](https://badges.gitter.im/macmanes-lab/Oyster_River_Protocol.svg)](https://gitter.im/macmanes-lab/Oyster_River_Protocol?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Official Repository of the Oyster River Protocol for Transcriptome Assembly

Please see http://oyster-river-protocol.readthedocs.io/en/latest/ and https://hackmd.io/s/SJhOQvkVm# for instructions about how to run and install. 

## Installation

See [INSTALL.md](INSTALL.md) for step-by-step install directions, covering both the one-command `make` installer and the manual, step-by-step alternative.

## Trinity read normalization

By default, oyster.py runs Trinity with `--no_normalize_reads`, i.e. read normalization is disabled. To let Trinity normalize reads instead, pass `--normalize-reads` on the command line:

    python3 oyster.py --read1 R1.fq.gz --read2 R2.fq.gz --mem 110 --cpu 24 --runout runname --strand RF --normalize-reads
