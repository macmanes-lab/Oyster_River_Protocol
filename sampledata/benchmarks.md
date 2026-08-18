# Benchmarks

## TIME2_SRR1789336_norm_py_parallel (ORP 3.0.0) -- 2026-08-17

Assembly: `/mnt/gpfs01/home/macmaneslab/macmanes/assemblies/TIME2_SRR1789336_norm_py_parallel.ORP.fasta`

Command:
```
oyster.py --normalize-reads --tpm-filt 1 --mem 670 --cpu 40 --read1 SRR1789336_1.fastq.gz --read2 SRR1789336_2.fastq.gz --max-parallel 2 --runout TIME2_SRR1789336_norm_py_parallel
```

| Metric | Value |
| --- | --- |
| BUSCO | C:88.0%[S:73.6%,D:14.4%],F:7.2%,M:4.8%,n:125 |
| TransRate score | 0.51827 |
| TransRate optimal score | 0.53225 |
| Unique genes (ORP) | 13772 |
| Unique genes (Trinity) | 13231 |
| Unique genes (SPAdes55) | 13493 |
| Unique genes (SPAdes75) | 12809 |
| Unique genes (Trans-ABySS) | 12632 |
| Reads mapped as proper pairs | 96.06% |

Step timing:
```
run_trimmomatic  00:04:21
run_rcorrector   00:14:25
run_transabyss   04:46:42
run_spades75     00:07:48
run_spades55     00:09:24
run_trinity      37:11:49
run_filtershort  00:00:07
run_orthofuser   00:05:15
orthofuser_branch 00:06:08
orthotransrate   00:26:36
merge_branch     00:26:37
makeorthout      00:00:49
orthofusing      00:05:35
diamond_orthomerged 00:00:08
diamond_transabyss 00:00:07
diamond_spades75 00:00:06
diamond_spades55 00:00:06
diamond_trinity  00:00:08
make_list5       00:00:00
posthack         00:00:04
cdhit            00:00:34
orp_diamond      00:00:07
salmon_index     00:00:21
salmon           00:00:29
secondfilter     00:00:54
busco            00:05:34
strandeval       00:01:34
transrate        00:15:33
TOTAL            38:28:11
```

## TIME2_SRR1789336_norm_py_NOparallel (ORP 3.0.0) -- 2026-08-17

Assembly: `/mnt/gpfs01/home/macmaneslab/macmanes/assemblies/TIME2_SRR1789336_norm_py_NOparallel.ORP.fasta`

Command:
```
oyster.py --normalize-reads --tpm-filt 1 --mem 670 --cpu 40 --read1 SRR1789336_1.fastq.gz --read2 SRR1789336_2.fastq.gz --max-parallel 1 --runout TIME2_SRR1789336_norm_py_NOparallel
```

| Metric | Value |
| --- | --- |
| BUSCO | C:91.2%[S:68.8%,D:22.4%],F:4.8%,M:4.0%,n:125 |
| TransRate score | 0.51346 |
| TransRate optimal score | 0.53753 |
| Unique genes (ORP) | 13765 |
| Unique genes (Trinity) | 13228 |
| Unique genes (SPAdes55) | 13493 |
| Unique genes (SPAdes75) | 12809 |
| Unique genes (Trans-ABySS) | 12632 |
| Reads mapped as proper pairs | 97.13% |

Step timing:
```
run_trimmomatic  00:05:11
run_rcorrector   00:14:11
run_trinity      36:51:56
run_transabyss   04:43:36
run_spades75     00:07:09
run_spades55     00:07:31
run_filtershort  00:00:07
run_orthofuser   00:30:57
orthofuser_branch 00:31:56
orthotransrate   00:23:41
merge_branch     00:23:41
makeorthout      00:01:08
orthofusing      00:07:54
diamond_orthomerged 00:00:09
diamond_transabyss 00:00:07
diamond_spades75 00:00:06
diamond_spades55 00:00:07
diamond_trinity  00:00:08
make_list5       00:00:00
posthack         00:00:04
cdhit            00:00:39
orp_diamond      00:00:08
salmon_index     00:00:21
salmon           00:00:30
secondfilter     00:01:02
busco            00:05:47
strandeval       00:01:33
transrate        00:14:27
TOTAL            43:39:49
```

## `--max-parallel 2` vs `--max-parallel 1` comparison

| Metric | parallel (2) | NOparallel (1) |
| --- | --- | --- |
| BUSCO complete | 88.0% | 91.2% |
| BUSCO single-copy | 73.6% | 68.8% |
| BUSCO duplicated | 14.4% | 22.4% |
| BUSCO fragmented | 7.2% | 4.8% |
| BUSCO missing | 4.8% | 4.0% |
| TransRate score | 0.51827 | 0.51346 |
| TransRate optimal score | 0.53225 | 0.53753 |
| Unique genes (ORP) | 13772 | 13765 |
| Unique genes (Trinity) | 13231 | 13228 |
| Unique genes (SPAdes55) | 13493 | 13493 |
| Unique genes (SPAdes75) | 12809 | 12809 |
| Unique genes (Trans-ABySS) | 12632 | 12632 |
| Reads mapped as proper pairs | 96.06% | 97.13% |
| Total wallclock time | 38:28:11 | 43:39:49 |

### Assembly metric deltas (NOparallel minus parallel)

| Metric | parallel (2) | NOparallel (1) | Delta |
| --- | --- | --- | --- |
| BUSCO complete | 88.0% | 91.2% | +3.2 pp |
| BUSCO single-copy | 73.6% | 68.8% | -4.8 pp |
| BUSCO duplicated | 14.4% | 22.4% | +8.0 pp |
| BUSCO fragmented | 7.2% | 4.8% | -2.4 pp |
| BUSCO missing | 4.8% | 4.0% | -0.8 pp |
| TransRate score | 0.51827 | 0.51346 | -0.00481 |
| TransRate optimal score | 0.53225 | 0.53753 | +0.00528 |
| Unique genes (ORP) | 13772 | 13765 | -7 |
| Unique genes (Trinity) | 13231 | 13228 | -3 |
| Unique genes (SPAdes55) | 13493 | 13493 | 0 |
| Unique genes (SPAdes75) | 12809 | 12809 | 0 |
| Unique genes (Trans-ABySS) | 12632 | 12632 | 0 |
| Reads mapped as proper pairs | 96.06% | 97.13% | +1.07 pp |

SPAdes55/75 and Trans-ABySS unique-gene counts are identical between runs (deterministic, single-threaded-per-step assemblers), so all of the ORP/Trinity/BUSCO/TransRate movement traces back to Trinity, whose result varies with concurrency (non-deterministic multi-threaded assembly). The BUSCO shift is the most notable: NOparallel gained 3.2pp complete but nearly all of that came from a jump in duplicated BUSCOs (+8.0pp) rather than single-copy (-4.8pp), which is a mixed signal rather than a clear quality win. `--max-parallel 2` finished 5h11m38s faster overall, mostly from `orthofuser_branch`/`merge_branch` and `run_orthofuser` running concurrently instead of sequentially (each roughly 5x longer under `--max-parallel 1`).

## test37 (`--max-parallel 2`) vs test36 (`--max-parallel 1`) -- 2026-08-18

Sanity-check run on the small CI sample dataset (`sampledata/test.1.fq.gz` / `test.2.fq.gz`) rather than a full-size assembly. No BUSCO/TransRate/unique-gene metrics were captured for this pair -- step timing only.

Command (test37, parallel):
```
oyster.py --strand RF --tpm-filt 0.2 --mem 15 --cpu 8 --read1 test.1.fq.gz --read2 test.2.fq.gz --max-parallel 2 --normalize-reads --runout test37
```

Command (test36, NOparallel):
```
oyster.py --strand RF --tpm-filt 0.2 --mem 15 --cpu 8 --read1 test.1.fq.gz --read2 test.2.fq.gz --max-parallel 1 --normalize-reads --runout test36
```

| Step | test37 (parallel 2) | test36 (parallel 1) |
| --- | --- | --- |
| run_trimmomatic | 00:00:01 | 00:00:01 |
| run_rcorrector | 00:00:01 | 00:00:01 |
| run_trinity_phase1 | 00:00:06 | 00:00:07 |
| run_transabyss | 00:00:08 | 00:00:08 |
| diamond_transabyss | 00:00:12 | 00:00:12 |
| run_spades55 | 00:00:01 | 00:00:01 |
| diamond_spades55 | 00:00:12 | 00:00:12 |
| run_spades75 | 00:00:01 | 00:00:01 |
| diamond_spades75 | 00:00:12 | 00:00:12 |
| run_trinity_phase2 | 00:00:58 | 00:01:00 |
| run_filtershort | 00:00:02 | 00:00:02 |
| run_orthofuser | 00:00:05 | 00:00:05 |
| orthofuser_branch | 00:00:05 | 00:00:05 |
| orthotransrate | 00:00:06 | 00:00:05 |
| merge_branch | 00:00:06 | 00:00:05 |
| makeorthout | 00:00:00 | 00:00:00 |
| orthofusing | 00:00:00 | 00:00:00 |
| diamond_orthomerged | 00:00:07 | 00:00:07 |
| diamond_trinity | 00:00:07 | 00:00:07 |
| make_list5 | 00:00:00 | 00:00:00 |
| posthack | 00:00:00 | 00:00:00 |
| cdhit | 00:00:00 | 00:00:00 |
| orp_diamond | 00:00:07 | 00:00:07 |
| salmon_index | 00:00:01 | 00:00:01 |
| salmon | 00:00:00 | 00:00:00 |
| secondfilter | 00:00:00 | 00:00:00 |
| busco | 00:01:17 | 00:01:19 |
| transrate | 00:00:05 | 00:00:04 |
| strandeval | 00:00:09 | 00:00:09 |
| TOTAL | 00:04:00 | 00:04:13 |

`--max-parallel 2` finished 13s faster. The step print order confirms the concurrency logic is engaging correctly: under `--max-parallel 2`, `run_orthofuser`/`orthofuser_branch` and `orthotransrate`/`merge_branch` interleave in the log (concurrent completion), while under `--max-parallel 1` they print strictly branch-by-branch (`orthotransrate`, `merge_branch`, then `run_orthofuser`, `orthofuser_branch`). Unlike the SRR1789336 benchmark above, each branch here only takes ~5s on this tiny sample dataset, so overlapping them saves only a few seconds rather than hours -- this run validates the scheduling behavior, not the wall-clock benefit at scale.
