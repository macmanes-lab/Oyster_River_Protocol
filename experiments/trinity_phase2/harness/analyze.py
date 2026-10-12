#!/usr/bin/env python3
"""Aggregate a replay run: per-stage CPU, per-job totals, java/python calls."""
import sys, os, re, glob, collections

work = sys.argv[1]
def cat(label):
    t = label.split()
    if not t: return "?"
    if t[0] == "bash" and "bowtie2" in label: return "bowtie2|samtools"
    b = os.path.basename(t[0])
    if b == "ParaFly":
        return "ParaFly(butterfly)" if "butterfly" in label.lower() else "ParaFly(quantifygraph)"
    if b == "inchworm":
        m = re.search(r"-K (\d+)", label)
        return "inchworm(super-read K=%s)" % m.group(1) if m and m.group(1) != "25" else "inchworm(main K=25)"
    if b in ("touch", "mkdir", "ln", "mv", "rm", "cat", "sed", "sort", "cp"): return "shell:" + b
    return b

stage = collections.defaultdict(lambda: [0, 0.0, 0.0, 0.0])  # n, wall, user, sys
jobs = []
for f in glob.glob(os.path.join(work, "prof", "*.tsv")):
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if p[1] == "CMD":
            s = stage[cat(p[5])]
            s[0] += 1; s[1] += float(p[2]); s[2] += float(p[3]); s[3] += float(p[4])
        elif p[1] == "TOTAL":
            s = stage["perl driver (self)"]
            s[0] += 1; s[2] += float(p[3]); s[3] += float(p[4])
for f in glob.glob(os.path.join(work, "time", "*.txt")):
    try:
        w, u, s_, m = open(f).read().split()[-4:]
        jobs.append((float(w), float(u) + float(s_), int(m), os.path.basename(f)))
    except Exception:
        pass
totcpu = sum(j[1] for j in jobs)
attributed = sum(v[2] + v[3] for v in stage.values())
print("jobs=%d  total job CPU=%.0f s  attributed=%.0f s  unattributed(backticks, perl startup)=%.0f s"
      % (len(jobs), totcpu, attributed, totcpu - attributed))
print("%-32s %7s %9s %9s %6s" % ("stage", "calls", "wall_s", "cpu_s", "cpu%"))
for k, v in sorted(stage.items(), key=lambda kv: -(kv[1][2] + kv[1][3])):
    c = v[2] + v[3]
    print("%-32s %7d %9.0f %9.0f %5.1f%%" % (k, v[0], v[1], c, 100 * c / totcpu if totcpu else 0))
# job size buckets
jobs.sort()
n = len(jobs)
if n:
    print("job wall quantiles:", " ".join("q%g=%.1f" % (q, jobs[min(n - 1, int(q * n))][0]) for q in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)))
    bycpu = sorted(jobs, key=lambda j: j[1])
    tiny = [j for j in jobs if j[0] < 15]
    print("jobs <15 s wall: %d (%.0f%%) holding %.0f%% of CPU" % (len(tiny), 100 * len(tiny) / n, 100 * sum(j[1] for j in tiny) / totcpu))
    print("top 2%% jobs by CPU hold %.0f%% of CPU" % (100 * sum(j[1] for j in bycpu[-max(1, n // 50):]) / totcpu))
jp = os.path.join(work, "java_py.tsv")
if os.path.exists(jp):
    J = collections.defaultdict(list)
    for line in open(jp):
        p = line.rstrip("\n").split("\t")
        if p[0] in ("JAVA", "PY"):
            try:
                J[p[0]].append((float(p[1]), float(p[2]) + float(p[3]), int(p[4])))
            except ValueError:
                pass
    for k, v in J.items():
        v.sort(key=lambda x: x[1])
        m = len(v)
        print("%s calls=%d wall=%.0f cpu=%.0f median_cpu=%.2f median_wall=%.2f median_rss_MB=%.0f cpu/wall=%.2f"
              % (k, m, sum(x[0] for x in v), sum(x[1] for x in v), v[m // 2][1],
                 sorted(x[0] for x in v)[m // 2], sorted(x[2] for x in v)[m // 2] / 1024,
                 sum(x[1] for x in v) / max(1e-9, sum(x[0] for x in v))))
        small = [x for x in v if x[0] < 2]
        print("   %s calls under 2 s wall: %d, their CPU=%.0f s" % (k, len(small), sum(x[1] for x in small)))
