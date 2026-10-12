#!/bin/bash
E=~/orp_envs/orp_trinity; export PATH=$E/bin:$PATH
cd ~/trinity_eng
for n in c1122:Fb_0/CBin_11 c3226:Fb_0/CBin_32; do c=${n%%:*}; p=${n#*:}
for tr in perf; do for jar in shipped unmod guard; do for rep in 1; do
  T=~/trinity_eng/src_detfix/Trinity; [ $tr = perf ] && T=~/trinity_eng/src/Trinity
  J=~/trinity_eng/src/Butterfly/Butterfly.jar; [ $jar = guard ] && J=~/trinity_eng/Butterfly.devel_guard.jar; [ $jar = unmod ] && J=~/trinity_eng/Butterfly.devel_unmod.jar
  d=/tmp/iso.$$.$c.$tr.$jar.$rep; mkdir -p $d; cp ~/assemblies/TIME2_SRR1789336_norm_py_5050parallel.trinity/read_partitions/$p/$c.trinity.reads.fa $d/
  (cd $d; CHRYSALIS_KMER_TABLE=$([ $tr = detfix ] && echo dense || echo auto) $T --single $c.trinity.reads.fa --output $c.trinity.reads.fa.out --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --no_version_check --bypass_java_version_check --no_distributed_trinity_exec --no_salmon --bfly_jar=$J >/dev/null 2>&1)
  echo "$c trinity=$tr jar=$jar rep=$rep $(grep -v ">" $d/$c.trinity.reads.fa.out.Trinity.fasta | md5sum | cut -c1-8)"; rm -rf $d
done; done; done; done
