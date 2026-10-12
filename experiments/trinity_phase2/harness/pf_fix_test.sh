#!/bin/bash
# test ParaFly configure under bioconda-like and plain flags, before/after fix
set -u
test_build() { # label, then env assignments
  local label=$1; shift
  make distclean >/dev/null 2>&1
  env "$@" ./configure --prefix=$(pwd) > cfg.log 2>&1
  local am=$(grep -E "^AM_CXXFLAGS" src/Makefile)
  env "$@" make > mk.log 2>&1
  echo "  $label | GOMP=$(nm -D src/ParaFly 2>/dev/null | grep -c GOMP_) | $am | $(grep 'option to support OpenMP' cfg.log) | errs=$(grep -c 'command not found' cfg.log)"
}
cd /tmp && rm -rf pf_fix && cp -r ~/trinity_eng/src_upstream/trinity-plugins/ParaFly pf_fix && cd pf_fix
echo "BEFORE (shipped configure):"
test_build plain
test_build bioconda LDFLAGS=-fopenmp CXXFLAGS=-O3
test_build trinity-makefile CFLAGS=-fopenmp CXXFLAGS=-fopenmp
python3 - <<'PY'
s = open("configure.ac").read()
old = "AC_OPENMP\nAC_SUBST([AM_CXXFLAGS], [-m64 $OPENMP_CXXFLAGS])"
assert old in s
s = s.replace(old, '''# AC_OPENMP compiles its probe with LDFLAGS, so -fopenmp there (as conda
# builds set it) reads as "none needed" while the real compiles, which do not
# see LDFLAGS, get no OpenMP: ParaFly then runs its commands one at a time.
save_LDFLAGS="$LDFLAGS"
LDFLAGS=""
AC_OPENMP
LDFLAGS="$save_LDFLAGS"
AC_SUBST([AM_CXXFLAGS], ["-m64 $OPENMP_CXXFLAGS"])''')
open("configure.ac", "w").write(s)
PY
autoreconf -fi > ar.log 2>&1 || { echo autoreconf failed; tail ar.log; }
echo "AFTER (fixed configure.ac, autoreconf):"
test_build plain
test_build bioconda LDFLAGS=-fopenmp CXXFLAGS=-O3
test_build trinity-makefile CFLAGS=-fopenmp CXXFLAGS=-fopenmp
