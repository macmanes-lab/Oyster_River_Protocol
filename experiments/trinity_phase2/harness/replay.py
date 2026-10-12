#!/usr/bin/env python3
"""Replay a seeded sample of real Trinity Phase 2 commands.

replay.py setup --ref REFDIR --n N --seed S --out SAMPLE.tsv
replay.py stage --sample SAMPLE.tsv --work WORK --trinity TRINITY [--extra "..."]
replay.py compare --sample SAMPLE.tsv --work WORK
"""
import argparse, os, random, re, shlex, shutil, sys, hashlib

def setup(a):
    lines = open(os.path.join(a.ref, "recursive_trinity.cmds")).read().splitlines()
    rnd = random.Random(a.seed)
    pick = sorted(rnd.sample(range(len(lines)), a.n))
    with open(a.out, "w") as fh:
        for i in pick:
            argv = shlex.split(lines[i])
            reads = argv[argv.index("--single") + 1]
            flags = [x for x in argv[1:]]
            # drop --single X --output Y, keep the rest verbatim
            j = flags.index("--single"); del flags[j:j + 2]
            j = flags.index("--output"); del flags[j:j + 2]
            rel = os.path.relpath(os.path.realpath(reads), os.path.realpath(a.ref))
            fh.write("%d\t%s\t%s\n" % (i, rel, " ".join(flags)))

def read_sample(path):
    for line in open(path):
        i, rel, flags = line.rstrip("\n").split("\t")
        yield int(i), rel, flags

def stage(a):
    os.makedirs(os.path.join(a.work, "prof"), exist_ok=True)
    os.makedirs(os.path.join(a.work, "time"), exist_ok=True)
    cmds = open(os.path.join(a.work, "cmds"), "w")
    for i, rel, flags in read_sample(a.sample):
        src = os.path.join(a.ref, rel)
        dst = os.path.join(a.work, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if not os.path.exists(dst):
            shutil.copyfile(src, dst)
        tag = "j%05d" % i
        if a.extra_strip:
            for f in a.extra_strip.split():
                flags = flags.replace(f, "")
        cmd = ("cd %s && TRINITY_PROF_LOG=%s /usr/bin/time -f '%%e\\t%%U\\t%%S\\t%%M' -o %s %s --single %s --output %s %s %s > %s.log 2>&1"
               % (shlex.quote(os.path.dirname(dst)),
                  os.path.join(a.work, "prof", tag + ".tsv"),
                  os.path.join(a.work, "time", tag + ".txt"),
                  a.trinity, shlex.quote(dst), shlex.quote(dst + ".out"),
                  flags, a.extra or "", shlex.quote(dst)))
        cmds.write(cmd + "\n")
    cmds.close()

def seqs(path):
    out = {}
    if not os.path.exists(path):
        return None
    name = None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]; out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}

def compare(a):
    same = diff = missing_both = 0
    bad = []
    for i, rel, flags in read_sample(a.sample):
        ref = seqs(os.path.join(a.ref, rel + ".out.Trinity.fasta"))
        new = seqs(os.path.join(a.work, rel + ".out.Trinity.fasta"))
        if ref is None and new is None:
            missing_both += 1
        elif ref == new:
            same += 1
        elif ref is not None and new is not None and sorted(ref.values()) == sorted(new.values()):
            same += 1
        else:
            diff += 1; bad.append(rel)
    print("identical=%d differ=%d no_output_either=%d" % (same, diff, missing_both))
    for b in bad[:20]:
        print("  DIFF", b)

p = argparse.ArgumentParser()
p.add_argument("mode")
p.add_argument("--ref", default=os.path.expanduser("~/assemblies/TIME2_SRR1789336_norm_py_5050parallel.trinity"))
p.add_argument("--n", type=int, default=1500)
p.add_argument("--seed", type=int, default=1)
p.add_argument("--out"); p.add_argument("--sample"); p.add_argument("--work")
p.add_argument("--trinity"); p.add_argument("--extra", default="")
p.add_argument("--extra-strip", default="")
a = p.parse_args()
{"setup": setup, "stage": stage, "compare": compare}[a.mode](a)
