#!/bin/bash
# find debugMes call sites whose guarding changes c1122's output
E=~/orp_envs/orp_trinity; export PATH=$E/bin:$PATH
B=~/trinity_eng/src/Butterfly/Butterfly
W=/tmp/bis.$$; mkdir -p $W; cd $W
cp ~/assemblies/TIME2_SRR1789336_norm_py_5050parallel.trinity/read_partitions/Fb_0/CBin_11/c1122.trinity.reads.fa .
~/trinity_eng/src/Trinity --single c1122.trinity.reads.fa --output c1122.trinity.reads.fa.out --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --no_cleanup --no_version_check --bypass_java_version_check --no_distributed_trinity_exec --no_salmon >/dev/null 2>&1
cd c1122.trinity.reads.fa.out
CMDS=chrysalis/butterfly_commands
out_md5() { # jar
  local s=""
  while read C; do CJ=$(echo "$C" | sed "s#-jar [^ ]*#-jar $1#"); bash -c "$CJ" >/dev/null 2>&1; G=$(echo "$C" | grep -o -- "-C [^ ]*" | cut -d" " -f2); s="$s$(md5sum < $G.allProbPaths.fasta)"; done < $CMDS
  echo "$s" | md5sum | cut -c1-8
}
REF=$(out_md5 ~/trinity_eng/Butterfly.devel_unmod.jar)
echo "reference (unmod) $REF"
cd $B; git checkout -q origin/devel 2>/dev/null
build() { # lo hi exclude -> jar path
  git -C $B checkout -q origin/devel -- src/src/TransAssembly_allProbPaths.java
  python3 ~/trinity_eng/guard_debug2.py $B/src/src/TransAssembly_allProbPaths.java $1 $2 "$3" > /dev/null
  ~/trinity_eng/build_bfly.sh $W/t.jar > /dev/null 2>&1
  git -C $B checkout -q origin/devel -- src/src/TransAssembly_allProbPaths.java
  echo $W/t.jar
}
test_range() { local j=$(build $1 $2 "$3"); (cd $W/c1122.trinity.reads.fa.out && out_md5 $j); }
EXCL=""
for round in 1 2 3 4 5; do
  all=$(test_range 0 100000 "$EXCL")
  echo "round $round: guard all except [$EXCL] -> $all"
  [ "$all" = "$REF" ] && { echo "NEUTRAL with exclusions: $EXCL"; break; }
  lo=0; hi=800   # smallest hi with [0,hi) changing output
  while [ $((hi - lo)) -gt 1 ]; do
    mid=$(( (lo + hi) / 2 ))
    r=$(test_range 0 $mid "$EXCL")
    if [ "$r" = "$REF" ]; then lo=$mid; else hi=$mid; fi
  done
  echo "culprit call index: $lo"
  EXCL=$(echo "$EXCL,$lo" | sed "s/^,//")
done
git -C $B checkout -q devel-guard
rm -rf $W
