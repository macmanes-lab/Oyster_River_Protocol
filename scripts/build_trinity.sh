#!/bin/bash
# Build Trinity from source into the orp_trinity env, with OpenMP.
#
#   conda run -n orp_trinity scripts/build_trinity.sh [TRINITY_DIR]
#
# Run it inside the orp_trinity env: the env's compilers, cmake, jellyfish
# and so on are what the build uses, and the env's bin/ is where Trinity gets
# linked. TRINITY_DIR defaults to software/trinityrnaseq beside this
# repository.
#
# Why this exists rather than `trinity=2.15.2` from bioconda: builds _4 to _6
# (May 2025 on) ship a ParaFly compiled without OpenMP, so Trinity's phase 2
# -- the per-component assemblies, ~93% of a Trinity run -- runs one component
# at a time whatever --CPU says (34h27m against 2h03m with OpenMP on
# SRR1789336). This builds upstream's release unmodified; the OpenMP check at
# the end is what stops a serial ParaFly from slipping through again.
#
# Re-running is safe: it fetches only if the pinned commit is missing,
# rebuilds only what changed, and relinks.
set -euo pipefail

# Keep in step with the pin in INSTALL.md. A commit rather than a branch, so
# an install is reproducible. Its Chrysalis, Inchworm and Butterfly
# submodules are pinned by the commit itself.
TRINITY_REPO=${TRINITY_REPO:-https://github.com/trinityrnaseq/trinityrnaseq.git}
TRINITY_COMMIT=${TRINITY_COMMIT:-633f3cb6d7e1764fab99cb4f5675bad40ff87f08}  # Trinity-v2.15.2

here=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
dest=${1:-$here/software/trinityrnaseq}

[ -n "${CONDA_PREFIX:-}" ] || {
    echo "run inside the orp_trinity env: conda run -n orp_trinity $0" >&2
    exit 1
}
for tool in git make cmake g++ jellyfish bowtie2 samtools salmon java perl python3; do
    command -v "$tool" >/dev/null || { echo "$tool not found in $CONDA_PREFIX" >&2; exit 1; }
done

if [ ! -d "$dest/.git" ]; then
    mkdir -p "$(dirname "$dest")"
    git clone "$TRINITY_REPO" "$dest"
fi
cd "$dest"
git cat-file -e "$TRINITY_COMMIT^{commit}" 2>/dev/null || git fetch origin
git checkout -q "$TRINITY_COMMIT"
git submodule sync -q
# Not --recursive: Butterfly's own submodule is a private repo and fails the clone.
git submodule update --init -q

# no_bamsifter: bamsifter is only used by genome-guided runs, which ORP never
# does, and its build needs autoconf's autoheader.
jobs=$(nproc 2>/dev/null || echo 4)
make -j "$jobs" no_bamsifter

# A serial ParaFly is the failure this build exists to prevent, and it is
# silent: the build succeeds and phase 2 just takes a day longer.
parafly=$dest/trinity-plugins/BIN/ParaFly
[ "$(nm -D "$parafly" | grep -c GOMP_)" -gt 0 ] || {
    echo "$parafly was built without OpenMP" >&2
    exit 1
}

# Trinity locates its PerlLib, Chrysalis, Butterfly and plugins from its own
# resolved path, so a symlink is enough to put it on the env's PATH -- and
# `conda run -n orp_trinity Trinity` keeps working unchanged.
ln -sf "$dest/Trinity" "$CONDA_PREFIX/bin/Trinity"

echo "Trinity $(git describe --always --tags) built at $dest"
echo "ParaFly OpenMP symbols: $(nm -D "$parafly" | grep -c GOMP_)"
echo "linked $CONDA_PREFIX/bin/Trinity"
