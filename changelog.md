### CHANGELOG

ORP Version 3.0.0 <- 2.4.0

- replace oyster.mk (GNU Make) with oyster.py (Python 3.6+), the new primary pipeline entrypoint; every 2.4.0 fix and behavior below (software versions, BUSCO OrthoDB v12.2, NORMALIZE_READS, TPM/cd-hit/diamond-prerequisite fixes, etc.) carries forward unchanged except where noted here
- run independent stages concurrently, capped by the new `--max-parallel` flag (default 2, splitting `--cpu`/`--mem` across however many run at once): the 4 assemblers (Trinity/rnaSPAdes55/rnaSPAdes75/Trans-ABySS); the orthofuser vs. merge/orthotransrate branches; transrate vs. strandeval
- upgrade OrthoFinder 2.5.2 -> 3.1.5
- update software versions: Trinity 2.15.2, SPAdes 4.3.0, BUSCO 6.1.0, DIAMOND 2.2.5, salmon 2.5.1, samtools 1.24, bwa 0.7.19, rcorrector 1.0.7, cd-hit 4.8.1, trimmomatic 0.41, seqtk 1.5, mcl 22.282
- migrate BUSCO lineage database from OrthoDB v10 (eukaryota_odb10) to OrthoDB v12.2 (eukaryota_odb12.2), fetched via BUSCO's own `--download` instead of a hardcoded tarball URL
- update Python to 3.14 in orp_env.yml
- add [INSTALL.md](INSTALL.md): step-by-step install directions covering both the one-command `make` installer and the manual, environment-by-environment alternative, plus how to `mamba env update` an existing `orp` env and the known Python-3.14-solve gotchas above
- run the 4 assemblers as two fixed resource lanes instead of splitting `--cpu`/`--mem` evenly across all of them: Trinity gets its own lane with `TRINITY_LANE_SHARE` (80%) for the whole run, since it dominates wall-clock by orders of magnitude and its thread count is fixed once it launches (no way to hand it more CPU after the other three finish). The remaining 20% runs a second lane -- rnaSPAdes55, rnaSPAdes75, and Trans-ABySS in sequence -- firing each assembly's diamond search immediately after it finishes rather than waiting for the orthofuser/merge stage, which those particular diamond searches never actually depended on. Previously the even split meant Trinity ran alone at half `--cpu` for the last ~36 of its ~37 hours once the other three finished, with the other half of the machine idle
- add NORMALIZE_READS option to oyster.mk, exposing Trinity's `--no_normalize_reads` flag on the command line (defaults to FALSE, preserving prior behavior)


Minor/bugfix

- pin `python=3.14` for `orp_spades` specifically -- fixes SPAdes reporting "Python version 3.6.8 is not supported," a known SPAdes bug (ablab/spades#1319) where it detects a stray system Python instead of its own env's interpreter. Initially pinned across every environment for consistency, but reverted everywhere else after `orp_transabyss`'s solve failed (no non-RC Python 3.14 build available for that dependency set): each conda env's Python is fully isolated from every other env, so there's no real benefit to forcing the newest version outside the one env with an actual bug, and doing so just adds solver risk for older/less-actively-maintained packages
- add file-mtime-based resumability: a step whose output already exists and is newer than its inputs is skipped, so a failed/interrupted run can be re-invoked without recomputing finished work
- deliberately keep diamond, orp_diamond, salmon, and BUSCO sequential at full `--cpu` -- these are CPU-bound tools where splitting cores to enable overlap measured out to no net benefit (diamond) or a real regression (BUSCO, whose ~2 minutes dwarfs the ~13s of transrate+strandeval it would otherwise share cores with)
- add `--version` and full `argparse`-based `--help`
- retry a failed step up to twice (60s apart) before giving up, since a failure on a cluster is often transient (node preemption, filesystem hiccup) rather than a real bug; `run_spades()`/`run_transabyss()` also clean up their non-empty output dir before each retry, since both tools refuse to reuse one rather than resuming (Trinity is exempt from both, since it resumes from its own checkpoints in-place and its ~37h runtime means blindly retrying a deterministic failure could triple the wall time before giving up). Also log immediately when either of the two assembler-lane threads fails, rather than staying silent until the other lane -- often Trinity, tens of hours out -- separately finishes
- add `orp_orthofinder`/`orthofinder` to `check()`'s preflight tool verification, which never covered OrthoFinder before since it used to be a bare downloaded binary rather than a conda env
- remove the `orthofuser.py`-based version-check in the Makefile (`orthopath`/`orthufuserversion`, and the `orthofuser` target's curl/tar logic): dead code left over from ORP's pre-2021 install process, when a MacManes-lab-patched fork of OrthoFinder named `orthofuser.py` supplied nucleotide-input support that upstream OrthoFinder didn't have yet. That patch was retired in 2021 once upstream OrthoFinder added its own `-d` DNA-input flag, but the Makefile's check for it (a path that stopped existing at that point) was never cleaned up until now
- consolidate installed conda environments from 11 down to 5 (`orp`, `orp_spades`, `orp_trinity`, `orp_busco`, `orp_transabyss`) to make install faster and less failure-prone: `orp_rcorrector`, `orp_trimmomatic`, `orp_cdhit`, `orp_diamond`, and `orp_salmon` are folded into the base `orp` env, and `orp_sam` is removed entirely (its samtools was a version-for-version duplicate of `orp`'s, and its bwa/seqtk were never actually used from that env -- every bwa/seqtk call already ran inside `orp_trinity`). `orp_spades`/`orp_trinity`/`orp_busco`/`orp_transabyss` are left fully isolated since each has its own complex or fragile dependency tree -- `orp_trinity` in particular pins an old `salmon=1.10.3` for Trinity's internal isoform filtering, which must never resolve against the pipeline's own modern `salmon=2.5.1` now living in `orp`
- fix orthomerged diamond "unique gene count" rule listing spades55.diamond.txt as a prerequisite twice while omitting spades75.diamond.txt entirely
- fix readcheck's error path so a bad READ1/READ2 (missing file, too-short reads) actually reports its own error message when aborting the build, instead of failing on an unrelated "command not found"
- fix cd-hit-est's memory cap being hardcoded to 5000MB regardless of the run's requested MEM budget
- replace the list3->list5 diamond reconciliation step's shell loop (one grep per gene against the full diamond output) with a single hash-based Python pass, scripts/build_list5.py
- replace the per-orthogroup best-contig selection (one grep of contigs.csv per orthogroup via orthout.done) with a single hash-based Python pass, scripts/pick_best_contigs.py, retiring the orthout.done intermediate
- isolate all orp_* conda env creation from the user's/system's default conda channels
- pin a compatible salmon inside orp_trinity, since an unpinned resolve was picking a salmon 2.x whose Rust rewrite dropped flags Trinity's internal isoform filtering still calls with
- fix the TRANSRATE check to test the bundled binary directly instead of a stale proxy check
- remove dead RCORR/RCORRDIR/BUSCO/BUSCODIR/BUSCODB and TRINITY_KMER variables left over from earlier refactors
- blank out hardcoded /home/ubuntu paths in software/config.ini
- fix strandeval's `hist` histogram call: it was reading unfiltered raw data instead of the filtered numeric column, and separately losing its `-p '#'` argument when re-executed through conda run's internal shell
- fix the software-installed banner (SALMON/BUSCO/SPADES/etc.) printing between every timed step instead of once at the start, caused by per-step timing re-parsing the whole makefile
- fix the BUSCO score in qualreport picking up a line from BUSCO's JSON summary instead of its text summary, after newer BUSCO started writing both a .txt and .json short_summary per run

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
