# Working notes

Cross-machine scratchpad so a session on either machine can pick up where
the other left off. Keep entries short; newest on top. Delete/trim once
stale. The dated history through the 4.1.0 release prep is in
[docs/notes-archive.md](docs/notes-archive.md); code comments cite its
entries by date.

## 4.1.0 release checklist (2026-10-17)

`version.txt` reads 4.1.0, the changelog has one condensed 4.1.0 section,
and the README no longer says 5.0. Not tagged yet. Before tagging `4.1.0`:

- ~~**Fix:** chowder ignores a changed `--seed`/`--assembly-order` on resume.~~
  Fixed with the release review's resume findings: `ingest` now takes
  `assemblies/working/<run>.ingest.settings` (order, seed, and each input's
  size and mtime) as an input, as the assemblers take theirs.
- **Release review (2026-10-17), fixed:** cleanup glob vs dotted run names,
  chowder cross-input contig IDs, resume stamps, relative read paths,
  assembly_eval quoting, twotrack's 20-subject cap, strandeval threads
  (changelog Fixes; `tests/test_run_integrity.py`).
  `tests/test_twotrack_select.py` passed on Premise in the orp env (real
  blastn). Still unexercised: the new samtools sort line (no `samtools
  view`, `-m 512M`) -- the full tier's `strandeval_standalone` and every
  full run's strandeval step cover it.
- ~~**Fix:** `--inchworm_cpu 10` hard-coded.~~ Not a real bug: Trinity
  2.15.2 caps `--inchworm_cpu` at `--CPU` itself (`Trinity` l.1232), so
  Inchworm never oversubscribed. ORP now passes `min(cpu, 10)` anyway so the
  logged command is accurate (`TRINITY_INCHWORM_CPU`, release_check case
  `trinity_inchworm_cpu`).
- **Test:** one `make` against a real conda on Premise (conda module) and
  one in a fresh container. The dev32 Makefile was only tested with stubs.
- **Test:** `tests/release_check.sbatch` full tier on the final commit.
  The quick tier passes on the laptop (9 pass, 1 skip without conda).
- **Test:** `--min-cov auto` end to end on real data.
- **Test:** the dev29 phase 2 grow fix (`clear_unfinished_trinity_components`)
  on the cluster; the three dev26 runs it would have saved weren't re-run.
- Then tag `4.1.0` on that commit.

## Premise ops

- `ssh premise` is set up on the `mdm2000` Mac (`~/.ssh/config`: key
  `~/.ssh/id_ed25519_premise`, user `macmanes`, connection reuse for 10
  min); the key goes on Premise with `ssh-copy-id -i
  ~/.ssh/id_ed25519_premise.pub premise`. Premise resets ssh connections opened in quick succession
  (`kex_exchange_identification: Connection reset`); the reuse should
  avoid it, otherwise wait and retry.
- `/mnt/home` is a symlink to `/mnt/gpfs01/home`; Trinity writes physical
  paths, so compare with `pwd -P`.
- After an outage, check for ghost jobs: a job can show RUNNING with no
  processes left (`srun --jobid <id> --overlap ps`).
- Exclude node117 and node118 (`--exclude=node117,node118`); Premise's
  repo copy of `oyster_mpi.sbatch` may lack it (`~/oyster_mpi.sbatch` has
  it).
