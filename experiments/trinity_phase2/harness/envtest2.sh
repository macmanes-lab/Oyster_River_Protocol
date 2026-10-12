#!/bin/bash
E=~/orp_envs/orp_trinity; export PATH=$E/bin:$PATH
for comp in c1122:Fb_0/CBin_11 c3226:Fb_0/CBin_32; do c=${comp%%:*}; p=${comp#*:}
d=/tmp/et2.$$.$c; mkdir -p $d; cp ~/assemblies/TIME2_SRR1789336_norm_py_5050parallel.trinity/read_partitions/$p/$c.trinity.reads.fa $d/
cd $d; ~/trinity_eng/src/Trinity --single $c.trinity.reads.fa --output $c.trinity.reads.fa.out --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --no_cleanup --no_version_check --bypass_java_version_check --no_distributed_trinity_exec --no_salmon >/dev/null 2>&1
A=$d/$c.trinity.reads.fa.out
for loc in C en_US.UTF-8; do for mode in default 0 2 3; do for jar in stable stable_guard; do
  s=""
  while read C; do
    CJ=$(echo "$C" | sed "s#-jar [^ ]*#-jar $HOME/trinity_eng/Butterfly.devel_$jar.jar#")
    [ $mode != default ] && CJ=$(echo "$CJ" | sed "s#^java #java -XX:+UnlockExperimentalVMOptions -XX:hashCode=$mode #")
    (cd $A && LC_ALL=$loc bash -c "$CJ" >/dev/null 2>&1); G=$(echo "$C" | grep -o -- "-C [^ ]*" | cut -d" " -f2); s="$s$(md5sum < $G.allProbPaths.fasta)"
  done < $A/chrysalis/butterfly_commands
  echo "$c LC_ALL=$loc hashCode=$mode jar=$jar $(echo $s | md5sum | cut -c1-8)"
done; done; done
rm -rf $d; done
