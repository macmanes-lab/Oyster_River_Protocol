#!/bin/bash
# Put the phase-2 fixes into an installed (bioconda) Trinity, with a backup.
#
#   install_into_env.sh parafly TRINITY_HOME   # ParaFly with OpenMP only: assembly unchanged
#   install_into_env.sh full    TRINITY_HOME   # + Trinity/Chrysalis/Butterfly patches
#   install_into_env.sh restore TRINITY_HOME   # put the backed-up files back
#
# TRINITY_HOME is the directory holding the Trinity script, e.g.
# ~/orp_envs/orp_trinity/bin. Built pieces come from ENG (default ~/trinity_eng
# on Premise): the OpenMP ParaFly and sparse-table Chrysalis binaries from
# src/, the devel + reproducible-hash + debug-guard Butterfly jar, and the
# trinityrnaseq patches 0001-0004 (0005 only touches ParaFly's sources).
set -euo pipefail

# replace by rename, so a ParaFly or Chrysalis binary that a running job has
# open keeps its old inode instead of being overwritten under it
put() { cp "$1" "$2.new.$$" && chmod --reference="$2" "$2.new.$$" && mv -f "$2.new.$$" "$2"; }

MODE=${1:?usage: $0 parafly|full|restore TRINITY_HOME}
TH=$(cd "${2:?usage: $0 parafly|full|restore TRINITY_HOME}" && pwd)
ENG=${ENG:-$HOME/trinity_eng}
BACKUP=$TH/../trinity_phase2_backup.tar

[ -f "$TH/Trinity" ] || { echo "no Trinity script in $TH" >&2; exit 1; }

FILES=(
    Trinity
    PerlLib/COMMON.pm
    PerlLib/Pipeliner.pm
    util/support_scripts/scaffold_iworm_contigs.pl
    Analysis/SuperTranscripts
    Chrysalis/bin
    Butterfly/Butterfly.jar
    trinity-plugins/BIN/ParaFly
)

if [ "$MODE" = restore ]; then
    [ -f "$BACKUP" ] || { echo "no backup at $BACKUP" >&2; exit 1; }
    tar -xf "$BACKUP" -C "$TH"
    rm "$BACKUP"
    echo "restored $TH from $BACKUP"
    exit 0
fi

if [ -f "$BACKUP" ]; then
    echo "$BACKUP exists: already installed? run '$0 restore $TH' first" >&2
    exit 1
fi
(cd "$TH" && tar -cf "$BACKUP" "${FILES[@]}")
echo "backed up originals to $BACKUP"

parafly=$ENG/src/trinity-plugins/BIN/ParaFly
nm -D "$parafly" | grep -q GOMP_ || { echo "$parafly has no OpenMP" >&2; exit 1; }
put "$parafly" "$TH/trinity-plugins/BIN/ParaFly"
echo "ParaFly: OpenMP build installed"

if [ "$MODE" = full ]; then
    for p in "$ENG"/patches/trinityrnaseq/000[1-4]-*.patch; do
        patch -d "$TH" -p1 --forward --no-backup-if-mismatch -s < "$p"
        echo "applied $(basename "$p")"
    done
    for b in "$ENG"/src/Chrysalis/bin/*; do put "$b" "$TH/Chrysalis/bin/$(basename "$b")"; done
    echo "Chrysalis: sparse k-mer table build installed"
    put "$ENG/Butterfly.devel_stable_guard.jar" "$TH/Butterfly/Butterfly.jar"
    echo "Butterfly: devel + reproducible hash + debug guard installed"
    perl -I"$TH/PerlLib" -c "$TH/Trinity"
elif [ "$MODE" != parafly ]; then
    echo "unknown mode $MODE" >&2
    exit 1
fi
