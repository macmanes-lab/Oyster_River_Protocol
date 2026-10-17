# Working notes

Cross-machine scratchpad so a session on either machine can pick up where
the other left off. Keep entries short; newest on top. Delete/trim once
stale. The dated history through the 4.1.0 release prep is in
[docs/notes-archive.md](docs/notes-archive.md); code comments cite its
entries by date.

## 4.1.0 release checklist (2026-10-17)

`version.txt` reads 4.1.0, the changelog has one condensed 4.1.0 section,
and the README no longer says 5.0. Not tagged yet. Before tagging `4.1.0`:

- **Fix:** chowder ignores a changed `--seed`/`--assembly-order` on resume.
  `ingest` is skipped when `<run>.ingest.done` and the ingested copies
  exist, so a re-run with a new seed reuses the old order while the report
  prints the new seed. Compare `ingest.done`'s order with the requested one
  and re-ingest when they differ.
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
