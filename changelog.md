### CHANGELOG

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
