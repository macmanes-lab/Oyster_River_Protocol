# Trans-ABySS under MPI

Trans-ABySS is the long pole of every ORP run since Trinity's phase 2 went
parallel (NOTES.md 2026-10-09): 8h55m on DRR031870, 3h52m on SRR866209 and
3h10m on SRR1138704, all at 10 cores, with Trinity's lane idle for 1-6 hours.
Most of that time is single-threaded: `abyss-pe` only runs its unitig
assembly in parallel in Bloom-filter mode or under MPI (`np`), and
`oyster.py` asks for neither (NOTES.md 2026-08-19). `transabyss --mpi N`
sets `np`, and `orp_transabyss` ships `ABYSS-P` and `mpirun`.

`run.sbatch` runs the two on the same reads and the same 16 cores:

| task | arm | command |
|---|---|---|
| 1 | `threads` | what `run_transabyss` runs today: `--threads 16` |
| 2 | `mpi` | the same plus `--mpi 16` |

```bash
sbatch --export=ALL,R1=<run>/rcorr/<run>.TRIM_1P.cor.fq.gz,R2=<run>/rcorr/<run>.TRIM_2P.cor.fq.gz,NAME=DRR031870 \
    experiments/transabyss_mpi/run.sbatch
```

Each arm writes, under `ta_mpi_<NAME>/<arm>/`:

- `summary.tsv`: wall time, peak RSS, minutes spent under 1.5 cores, contig
  and base counts;
- `cpu.tsv`: cores used per minute and the binaries using them, which shows
  where the single-threaded stretch is and whether MPI removed it;
- `eval/reports/qualreport.<NAME>.<arm>`: pyTransRate, BUSCO and strandeval,
  from `assembly_eval.py`.

## What decides it

- **Time**: the `threads` arm also gives Trans-ABySS's time at 16 cores, which
  no ORP run has measured yet.
- **Output**: MPI is worth adopting only if the assemblies are close.
  `ABYSS-P` is a different implementation of the same assembler, so expect
  similar rather than identical output. Compare the two qualreports
  (`diff` works, the layout is the same) and the contig counts.

If `mpirun` refuses to start 16 ranks under a one-task allocation, the log
says so in the first minute. The script sets OpenMPI's oversubscribe options
for that, but a site build can override them; resubmitting with
`--ntasks=16 --cpus-per-task=1` is the other way round it.
