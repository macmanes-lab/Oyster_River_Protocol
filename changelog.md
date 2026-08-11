### CHANGELOG

ORP Version 2.4.0 <- 2.3.3

- update software versions: Trinity 2.15.2, SPAdes 4.3.0, BUSCO 6.1.0, DIAMOND 2.2.5, salmon 2.5.1, samtools 1.24, bwa 0.7.19, rcorrector 1.0.7, cd-hit 4.8.1, trimmomatic 0.41, seqtk 1.5, mcl 22.282
- migrate BUSCO lineage database from OrthoDB v10 (eukaryota_odb10) to OrthoDB v12.2 (eukaryota_odb12.2), fetched via BUSCO's own `--download` instead of a hardcoded tarball URL
- fix orthomerged diamond "unique gene count" rule listing spades55.diamond.txt as a prerequisite twice while omitting spades75.diamond.txt entirely
- fix readcheck's error path so a bad READ1/READ2 (missing file, too-short reads) actually reports its own error message when aborting the build, instead of failing on an unrelated "command not found"
- update Python to 3.14 in orp_env.yml
- fix cd-hit-est's memory cap being hardcoded to 5000MB regardless of the run's requested MEM budget
- replace the list3->list5 diamond reconciliation step's shell loop (one grep per gene against the full diamond output) with a single hash-based Python pass, scripts/build_list5.py
- replace the per-orthogroup best-contig selection (one grep of contigs.csv per orthogroup via orthout.done) with a single hash-based Python pass, scripts/pick_best_contigs.py, retiring the orthout.done intermediate

ORP Version 2.3.2 <- 2.3.1

- update docker base image to Ubuntu 20.04 and new docker image
- update software version for Orthofinder 2,5,2, BUSCO 5.1.2, Trinity 2.12, SPAdes 3.15.2, diamond 2.0.8, RCORRECTOR 1.04
- handle a corner case there TPM filter does not remove any transcripts, which in previous versions caused a crash.
-  change to use mamba during install

ORP Version 2.3.0 <- 2.2.8

- remove Python 2.x dependencies!
- Update BUSCO to 4.0.0 (and to OrthoDB v10 databases)
- Update Spades to 3.14
- Update Trinity to 2.8.5
- Update Orthofuser to be in line with OrthoFinder 2.3.9
- Update the report generation script to use BUSCO4
- Fix bug where crash if TPM_FILT flag was not used.
- Fix bug where 1st low expression transcript was being removed.
- several small bug fixes which were not likely exposed in real datasets.  
