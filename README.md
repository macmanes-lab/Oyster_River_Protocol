# Oyster River Protocol

[![Join the chat at https://gitter.im/macmanes-lab/Oyster_River_Protocol](https://badges.gitter.im/macmanes-lab/Oyster_River_Protocol.svg)](https://gitter.im/macmanes-lab/Oyster_River_Protocol?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Official Repository of the Oyster River Protocol for Transcriptome Assembly

Please see http://oyster-river-protocol.readthedocs.io/en/latest/ and https://hackmd.io/s/SJhOQvkVm# for instructions about how to run and install. 

## Trinity read normalization

By default, oyster.mk runs Trinity with `--no_normalize_reads`, i.e. read normalization is disabled. To let Trinity normalize reads instead, pass `NORMALIZE_READS=TRUE` on the command line:

    oyster.mk READ1= READ2= MEM=110 CPU=24 RUNOUT=runname STRAND=RF NORMALIZE_READS=TRUE
