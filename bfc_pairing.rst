============================
BFC and read pairing info
============================

There is what I consider a bug in bfc, read pairing information is stipped when using bfc to correct fastQ reads in the Illumina 1.8 format or newer. This has been described here: https://github.com/lh3/bfc/issues/8 and https://github.com/macmanes-lab/Oyster_River_Protocol/issues/2.

The bug requires, assuming you want to interleave *and* de-interleave reads, which I think most people will want to do, that you place the ``/1`` and ``/2`` tags in the fastQ headers. To do this, follow the below steps.

::

  sed 's_ _/1 _g' file.1.fq > edited_file.1.fq

and to for the /2 read file

::

  sed 's_ _/2 _g' file.2.fq > edited_file.2.fq
  
  Then supply the edited files to seqtk for merging, bfc for correction, and so on...
