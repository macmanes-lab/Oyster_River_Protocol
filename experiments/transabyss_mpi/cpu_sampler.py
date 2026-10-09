#!/usr/bin/env python3
"""cpu_sampler.py - how many cores a process tree uses, minute by minute.

Usage:
    cpu_sampler.py ROOT_PID OUT.tsv [INTERVAL_SECONDS]

Every INTERVAL (default 60) seconds, walks /proc for ROOT_PID and all its
descendants, and writes one line: the time, the cores the tree used over the
interval (CPU seconds / wall seconds), and the three commands that used the
most. Exits when ROOT_PID is gone.

Unlike ps's %CPU (a lifetime average) this is the use over each interval, so
a stretch where Trans-ABySS sits on one core shows up as ~1.0 against the
binary doing it (ABYSS, ABYSS-P, AdjList, ...). A process that starts and
exits inside one interval is missed; at 60 s that is noise for a run of hours.
Linux only.
"""

import os
import sys
import time

CLK_TCK = os.sysconf("SC_CLK_TCK")


def read_stat(pid):
    """(ppid, comm, starttime, cpu ticks) for pid, or None if it has gone."""
    try:
        with open("/proc/%d/stat" % pid) as f:
            data = f.read()
    except OSError:
        return None
    # comm is in parentheses and may itself contain spaces or ')'.
    lpar, rpar = data.index("("), data.rindex(")")
    comm = data[lpar + 1:rpar]
    fields = data[rpar + 2:].split()
    # fields[0] is field 3 (state) of proc(5): ppid is 4, utime 14, stime 15,
    # starttime 22.
    return int(fields[1]), comm, int(fields[19]), int(fields[11]) + int(fields[12])


def tree(root):
    """{(pid, starttime): (comm, ticks)} for root and every descendant."""
    stats = {}
    for name in os.listdir("/proc"):
        if name.isdigit():
            st = read_stat(int(name))
            if st:
                stats[int(name)] = st
    children = {}
    for pid, (ppid, _, _, _) in stats.items():
        children.setdefault(ppid, []).append(pid)
    out, todo = {}, [root]
    while todo:
        pid = todo.pop()
        if pid in stats:
            _, comm, start, ticks = stats[pid]
            out[(pid, start)] = (comm, ticks)
        todo.extend(children.get(pid, ()))
    return out


def main():
    root, path = int(sys.argv[1]), sys.argv[2]
    interval = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0
    prev, t_prev = tree(root), time.time()
    with open(path, "w") as out:
        out.write("time\tcores\ttop_commands\n")
        while os.path.exists("/proc/%d" % root):
            time.sleep(interval)
            cur, t_cur = tree(root), time.time()
            used = {}
            for key, (comm, ticks) in cur.items():
                # A process new this interval counts all its ticks.
                delta = ticks - prev.get(key, (comm, 0))[1]
                used[comm] = used.get(comm, 0) + delta
            wall = (t_cur - t_prev) * CLK_TCK
            total = sum(used.values()) / wall
            top = sorted(used.items(), key=lambda kv: -kv[1])[:3]
            top_s = " ".join("%s=%.1f" % (c, t / wall) for c, t in top if t > 0)
            out.write("%s\t%.2f\t%s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), total, top_s))
            out.flush()
            prev, t_prev = cur, t_cur


if __name__ == "__main__":
    main()
