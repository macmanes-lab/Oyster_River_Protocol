#!/bin/bash
E=~/orp_envs/orp_trinity; export PATH=$E/bin:$PATH
R=~/assemblies/TIME2_SRR1789336_norm_py_5050parallel.trinity
cd $R; find read_partitions -name "*.trinity.reads.fa" -printf "%s %p\n" | sort -rn | head -5 > /tmp/heavy.$$
cat /tmp/heavy.$$
JAR=~/trinity_eng/Butterfly.stablehash_guard.jar
while read sz rel; do
  n=$(basename $rel .trinity.reads.fa); d=/tmp/hv.$$.$n; mkdir -p $d; cp $R/$rel $d/
  ( cd $d; s=$(date +%s); ~/trinity_eng/src/Trinity --single $n.trinity.reads.fa --output $n.trinity.reads.fa.out --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --no_cleanup --no_version_check --bypass_java_version_check --no_distributed_trinity_exec --no_salmon --bfly_jar=$JAR >log 2>&1; echo "$n reads_bytes=$sz whole_job=$(( $(date +%s)-s ))s graphs=$(wc -l < $n.trinity.reads.fa.out/chrysalis/butterfly_commands)"
    cd $n.trinity.reads.fa.out
    for mode in C1 FULL; do
      t0=$(date +%s.%N); u=0
      while read C; do
        [ $mode = FULL ] && C=$(echo "$C" | sed "s/-XX:TieredStopAtLevel=1//")
        /usr/bin/time -a -o /tmp/hvt.$$.$mode -f "%e %U" bash -c "$C >/dev/null 2>&1"
      done < chrysalis/butterfly_commands
      awk -v n=$n -v m=$mode "{w+=\$1; u+=\$2} END{printf \"  %s %s butterfly total wall=%.1fs cpu=%.1fs\n\", n, m, w, u}" /tmp/hvt.$$.$mode; rm -f /tmp/hvt.$$.$mode
      G=$(find chrysalis/Component_bins -name "*allProbPaths.fasta" | sort | xargs cat | md5sum | cut -c1-8); echo "  $mode output $G"
    done )
  rm -rf $d
done < /tmp/heavy.$$
