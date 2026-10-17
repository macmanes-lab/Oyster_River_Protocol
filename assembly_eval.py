#!/usr/bin/env python3
"""assembly_eval.py - score any transcriptome assembly the way the ORP does.

Usage:
    assembly_eval.py --assembly my.fasta --read1 R1.fq.gz --read2 R2.fq.gz \\
                     --cpu 24 --mem 110 --runout myrun

Runs the three evaluations an ORP run ends with, on one assembly you
already have, and writes the same quality report:

  * pyTransRate  contig-level read-mapping score (PYTRANSRATE SCORE, OPTIMAL),
  * BUSCO        gene-set completeness (--lineage, default eukaryota_odb12.2),
  * strandeval   strand-specificity of the library against the assembly,
                 plus the proper-pair rate from the same mapping.

Nothing is assembled, merged or filtered, and the reads are used as given
(no trimming or error correction): pass the pair you want scored, normally
the one the assembly was built from. The three steps are independent, so they
run at the same time, --cpu and --mem split between them.

Output, under --outdir (default ./<runout>_eval):

    reports/qualreport.<run>              the report, in oyster.py's layout
    reports/pytransrate_<run>/            pyTransRate's output
    reports/run_<run>.ORP/                BUSCO's output
    reports/<run>.strandeval_summary.txt  the strand histogram
    reports/<run>.flagstat                samtools flagstat of the sampled mapping
    logs/<step>.log                       each tool's stdout/stderr

The report has no UNIQUE GENES lines: those come from ORP's per-assembler
diamond searches, and a lone assembly has no assemblers to compare. Everything
else is the line oyster.py writes, in the same order, so reports from the two
diff cleanly.

A step whose reports/<run>.<step>.done sentinel exists is skipped, so an
interrupted run resumes where it stopped; --force starts over. Tools run in
the same conda envs as oyster.py (orp, orp_busco, orp_trinity); --env-* flags
rename them, or --no-conda runs whatever is on PATH.
"""

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
import textwrap
import threading
import time
from pathlib import Path

from oyster import RED, RESET, line_buffer_stdio, strand_exam

HERE = Path(__file__).resolve().parent

STEPS = ("pytransrate", "busco", "strandeval")

#: Steps whose tool output is echoed to the screen as well as logged, each
#: line prefixed with the step so the two concurrent streams stay legible.
#: strandeval is a handful of short bwa/samtools calls whose chatter says
#: nothing a reader needs; it stays in its log.
LOUD_STEPS = ("pytransrate", "busco")

#: One lock for every echoed line, so two tools never interleave mid-line.
_PRINT_LOCK = threading.Lock()


def read_version():
    try:
        return (HERE / "version.txt").read_text().strip()
    except OSError:
        return "unknown"


class Eval:
    def __init__(self, args):
        self.args = args
        self.assembly = Path(args.assembly).resolve()
        self.read1 = Path(args.read1).resolve()
        self.read2 = Path(args.read2).resolve()
        self.runout = args.runout or self.assembly.name.split(".")[0]
        self.cpu = args.cpu
        self.mem = args.mem
        self.lineage = args.lineage
        self.outdir = Path(args.outdir).resolve() if args.outdir else Path.cwd() / f"{self.runout}_eval"
        self.reports = self.outdir / "reports"
        self.logs = self.outdir / "logs"
        self.work = self.outdir / "work"
        self.version = read_version()
        self.pytransrate_args = shlex.split(args.pytransrate_args or "")
        self.busco_config = Path(args.busco_config).resolve() if args.busco_config else (
            HERE / "software" / "config.ini")
        self.busco_dbs = Path(args.busco_download_path).resolve() if args.busco_download_path else (
            HERE / "busco_dbs")
        self.strand_script = Path(args.examine_strand).resolve() if args.examine_strand else (
            HERE / "scripts" / "examine_strand.pl")
        # Quoted, as in oyster.py, so the echoed line can be pasted back.
        self.run_cmd = " ".join(
            shlex.quote(a) for a in [Path(sys.argv[0]).name] + sys.argv[1:])

    def welcome(self):
        """Say what is being scored, with what, and where the output goes."""
        print(RED)
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("                    _.-~~~~~-._")
        print("                 .-~   o   o   ~-.        OYSTER RIVER PROTOCOL")
        print("                (   .-'~~~~~'-.   )         .--.")
        print("                 '-.___________.-'         ( () )--  assembly_eval")
        print("                     '-.___.-'              '--'")
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + RESET + "\n")
        print(f"    version   {self.version} (assembly_eval)")
        print(textwrap.fill(
            self.run_cmd, width=78, break_long_words=False, break_on_hyphens=False,
            initial_indent="    command   ", subsequent_indent=" " * 14))
        print()
        print(f"    assembly  {self.assembly}")
        print(f"    read1     {self.read1}")
        print(f"    read2     {self.read2}")
        print(f"    runout    {self.runout}")
        print(f"    outdir    {self.outdir}")
        print(f"    cpu/mem   {self.cpu} threads, {self.mem}G (split across the steps)")
        print(f"    steps     {', '.join(self.steps) or 'none'}"
              + (f"   (skipped: {', '.join(self.args.skip)})" if self.args.skip else ""))
        if "busco" in self.steps:
            print(f"    lineage   {self.lineage}")
        if self.pytransrate_args:
            print(f"    pytransrate args  {' '.join(self.pytransrate_args)}")
        print(f"    conda     {'off (tools from PATH)' if self.args.no_conda else f'envs {self.args.env}, {self.args.env_busco}, {self.args.env_trinity}'}")
        print("\n    This scores an assembly you already have, NOT a full ORP run: nothing")
        print("    is assembled, trimmed, corrected or filtered. The reads are used as")
        print("    given, and the report has no UNIQUE GENES lines.\n")

    # -- plumbing ------------------------------------------------------------

    def conda_cmd(self, env, cmd):
        cmd = [str(c) for c in cmd]
        if self.args.no_conda:
            return cmd
        return ["conda", "run", "--no-capture-output", "-n", env, *cmd]

    def run(self, step, cmd, cwd=None, stdout=None, env=None):
        """Run cmd, appending its output to logs/<step>.log (or `stdout`).

        For a LOUD_STEPS step, stdout and stderr are also echoed to the
        screen line by line, prefixed with `[step]`.
        """
        cmd = [str(c) for c in cmd]
        with open(self.logs / f"{step}.log", "a") as log:
            log.write(f"$ {' '.join(cmd)}\n")
            log.flush()
            if step not in LOUD_STEPS or stdout is not None:
                subprocess.run(cmd, cwd=cwd, env=env, check=True,
                               stdout=stdout or log, stderr=log)
                return
            with _PRINT_LOCK:
                print(f"[{step}] + {' '.join(cmd)}")
            proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, universal_newlines=True,
                                    errors="replace", bufsize=1)
            for line in proc.stdout:
                log.write(line)
                log.flush()
                with _PRINT_LOCK:
                    print(f"[{step}] {line.rstrip()}")
            if proc.wait():
                raise subprocess.CalledProcessError(proc.returncode, cmd)

    def sentinel(self, step):
        return self.reports / f"{self.runout}.{step}.done"

    def adopt_pre_rename_outputs(self):
        """Carry a run from before the step was renamed transrate -> pytransrate.

        Its output dir, sentinel and log are renamed in place, as oyster.py's
        adopt_pre_rename_reports does, so a resumed run doesn't score the
        assembly all over again.
        """
        for old, new in (
            (self.reports / f"transrate_{self.runout}", self.reports / f"pytransrate_{self.runout}"),
            (self.sentinel("transrate"), self.sentinel("pytransrate")),
            (self.logs / "transrate.log", self.logs / "pytransrate.log"),
        ):
            if old.exists() and not new.exists():
                old.rename(new)
                print(f"[resume] renamed {old} to {new}")

    def preflight(self):
        for label, path in (("assembly", self.assembly), ("--read1", self.read1),
                            ("--read2", self.read2)):
            if not path.is_file():
                sys.exit(f"{label} not found: {path}")
        if "strandeval" in self.steps and not self.strand_script.is_file():
            sys.exit(f"examine_strand.pl not found: {self.strand_script} (--examine-strand)")
        if self.args.no_conda:
            tools = {"pytransrate": ["pytransrate"], "busco": ["busco"],
                     "strandeval": ["bwa", "seqtk", "samtools", "perl", "hist"]}
            missing = [t for s in self.steps for t in tools[s] if shutil.which(t) is None]
            if missing:
                sys.exit(f"--no-conda: not on PATH: {', '.join(missing)}")
        elif shutil.which("conda") is None:
            sys.exit("conda not found on PATH (or pass --no-conda to use tools on PATH)")

    # -- steps ---------------------------------------------------------------

    def pytransrate(self, cpu, mem):
        outdir = self.reports / f"pytransrate_{self.runout}"
        if outdir.exists():
            shutil.rmtree(outdir)
        mem_args = []
        if not any(a.split("=")[0] in ("--max-memory", "--mem") for a in self.pytransrate_args):
            mem_args = ["--max-memory", f"{mem}G"]
        self.run("pytransrate", self.conda_cmd(self.args.env, [
            "pytransrate", "-o", outdir, "-a", self.assembly,
            "--left", self.read1, "--right", self.read2, "-t", cpu,
            *mem_args, *self.pytransrate_args]))

    def busco(self, cpu, mem):
        name = f"run_{self.runout}.ORP"
        # BUSCO refuses an existing -o directory, and writes it under cwd.
        target = self.reports / name
        if target.exists():
            shutil.rmtree(target)
        cmd = ["busco", "--offline", "--lineage", self.lineage,
               "-i", self.assembly, "-m", "transcriptome", "--cpu", cpu, "-o", name]
        if self.busco_dbs.is_dir():
            cmd += ["--download_path", self.busco_dbs]
        if self.busco_config.is_file():
            cmd += ["--config", self.busco_config]
        env = os.environ.copy()
        if self.busco_config.is_file():
            env["BUSCO_CONFIG_FILE"] = str(self.busco_config)
        self.run("busco", self.conda_cmd(self.args.env_busco, cmd),
                 cwd=self.reports, env=env)

    def trinity_perllib(self):
        res = subprocess.run(
            self.conda_cmd(self.args.env_trinity, ["which", "Trinity"]),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, check=True)
        return Path(res.stdout.strip()).resolve().parent / "PerlLib"

    def strandeval(self, cpu, mem):
        """oyster.py's strand_exam, so the histogram is the one a full run
        draws, with its paths quoted and its threads kept to `cpu`; only the
        envs (or PATH, under --no-conda) and the script are this tool's."""
        work = self.work / "strandeval"
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        envs = {"orp": self.args.env, "orp_trinity": self.args.env_trinity}
        summary = strand_exam(
            lambda cmd, cwd=None: self.run("strandeval", cmd, cwd=cwd),
            self.assembly, self.read1, self.read2, self.runout, work, cpu,
            flagstat_out=self.reports / f"{self.runout}.flagstat",
            keep_scratch=True,
            wrap=lambda env, cmd: self.conda_cmd(envs[env], cmd),
            strand_script=self.strand_script, perllib=self.trinity_perllib(),
        )
        (self.reports / f"{self.runout}.strandeval_summary.txt").write_text(summary)
        if not self.args.no_cleanup:
            shutil.rmtree(work, ignore_errors=True)

    # -- report --------------------------------------------------------------

    def reportgen(self):
        run = self.runout
        lines = []

        def emit(label, value):
            text = f"{label}      {value}"
            print(text)
            lines.append(text)

        header = f"*****  QUALITY REPORT FOR: {run} using assembly_eval version {self.version} ****"
        print(f"\n\n{header}")
        lines.append(header)
        print(f"\n*****  THE ASSEMBLY SCORED IS: {self.assembly} **** \n")

        busco_line = ""
        for p in (self.reports / f"run_{run}.ORP").rglob("short*.txt"):
            for line in open(p):
                if re.match(r"^\s*C:[0-9]", line):
                    busco_line = line.strip()
        emit("*****  BUSCO SCORE ~~~~~~~~~~~~~~~~~~~~~~>", busco_line)

        csv_path = next((self.reports / f"pytransrate_{run}").rglob("assemblies.csv"), None)
        rows = list(csv.reader(open(csv_path))) if csv_path else []
        score = rows[1][36] if len(rows) > 1 and len(rows[1]) > 36 else ""
        optimal = rows[1][37] if len(rows) > 1 and len(rows[1]) > 37 else ""
        emit("*****  PYTRANSRATE SCORE ~~~~~~~~~~~~~~~~>     ", score)
        emit("*****  PYTRANSRATE OPTIMAL SCORE ~~~~~~~~>     ", optimal)

        proper_pairs = ""
        flagstat = self.reports / f"{run}.flagstat"
        if flagstat.exists():
            for line in open(flagstat):
                if "properly paired" in line:
                    fields = line.split()
                    if len(fields) > 5:
                        proper_pairs = fields[5].lstrip("(")
        emit("*****  READS MAPPED AS PROPER PAIRS ~~~~~>     ", proper_pairs)
        print(" \n")

        summary_path = self.reports / f"{run}.strandeval_summary.txt"
        if summary_path.exists():
            summary = summary_path.read_text()
            print(summary)
            lines.append(summary.rstrip("\n"))

        (self.reports / f"qualreport.{run}").write_text("\n".join(lines) + "\n")

    # -- orchestration -------------------------------------------------------

    def main(self):
        start = time.time()
        self.steps = [s for s in STEPS if s not in self.args.skip]
        self.welcome()
        self.preflight()
        for d in (self.reports, self.logs, self.work):
            d.mkdir(parents=True, exist_ok=True)
        self.adopt_pre_rename_outputs()
        if self.args.force:
            for s in STEPS:
                # Not unlink(missing_ok=True): that is 3.8+, and this runs
                # under the cluster's system python3, 3.6.
                if self.sentinel(s).exists():
                    self.sentinel(s).unlink()

        todo = [s for s in self.steps if not self.sentinel(s).exists()]
        for s in self.steps:
            if s not in todo:
                print(f"[{s}] already done, skipping (--force to redo)")

        # The steps are independent; --cpu and --mem are split across the ones
        # that run so the box is not oversubscribed. strandeval's mapping is
        # light and BUSCO scales poorly past a few threads, so pytransrate,
        # the slow one, gets the larger share.
        weights = {"pytransrate": 2, "busco": 1, "strandeval": 1}
        total = sum(weights[s] for s in todo) or 1
        failures = []

        def worker(step):
            share = weights[step] / total
            cpu = max(1, int(self.cpu * share))
            mem = max(1, int(self.mem * share))
            t0 = time.time()
            print(f"[{step}] started (cpu {cpu}, mem {mem}G); log: {self.logs / (step + '.log')}")
            try:
                getattr(self, step)(cpu, mem)
                self.sentinel(step).touch()
                print(f"[{step}] finished in {int(time.time() - t0)}s")
            except (subprocess.CalledProcessError, OSError) as e:
                failures.append((step, e))
                print(f"[{step}] FAILED: {e}  (see {self.logs / (step + '.log')})", file=sys.stderr)

        threads = [threading.Thread(target=worker, args=(s,)) for s in todo]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        self.reportgen()
        print(f"\nQuality report: {self.reports / ('qualreport.' + self.runout)}")
        print(f"Wall time: {int(time.time() - start)}s")
        if failures:
            sys.exit(f"failed: {', '.join(s for s, _ in failures)} "
                     "(report written with those lines blank; fix and re-run to resume)")


def parse_args():
    p = argparse.ArgumentParser(
        description="Score a transcriptome assembly with pytransrate, BUSCO and strandeval.")
    p.add_argument("--assembly", required=True, help="assembly fasta to score (may not be gzipped)")
    p.add_argument("--read1", required=True, help="left reads (fastq, optionally gzipped)")
    p.add_argument("--read2", required=True, help="right reads")
    p.add_argument("--runout", help="run name used in file names (default: assembly filename stem)")
    p.add_argument("--outdir", help="output directory (default: ./<runout>_eval)")
    p.add_argument("--cpu", type=int, default=os.cpu_count() or 8, help="total threads (default: all)")
    p.add_argument("--mem", type=int, default=32, help="total memory in GB (default: 32)")
    p.add_argument("--lineage", default="eukaryota_odb12.2",
                   help="BUSCO lineage (default: eukaryota_odb12.2)")
    # "transrate" is the step's name before it was renamed, still accepted.
    p.add_argument("--skip", nargs="*", default=[], choices=STEPS + ("transrate",), metavar="STEP",
                   help=f"steps not to run: {', '.join(STEPS)}")
    p.add_argument("--force", action="store_true", help="redo steps already marked done")
    p.add_argument("--no-cleanup", action="store_true", help="keep strandeval's BAM and bwa index")
    p.add_argument("--pytransrate-args", default="",
                   help="extra pytransrate arguments, as one quoted string")
    p.add_argument("--busco-config", help="BUSCO config.ini (default: software/config.ini if present)")
    p.add_argument("--busco-download-path",
                   help="BUSCO lineage dataset directory (default: busco_dbs/ beside this script if present)")
    p.add_argument("--examine-strand",
                   help="path to examine_strand.pl (default: scripts/ beside this script)")
    p.add_argument("--no-conda", action="store_true", help="run tools from PATH, not conda envs")
    p.add_argument("--env", default="orp", help="conda env with samtools (default: orp)")
    p.add_argument("--env-busco", default="orp_busco", help="conda env with busco (default: orp_busco)")
    p.add_argument("--env-trinity", default="orp_trinity",
                   help="conda env with Trinity, bwa, seqtk, hist (default: orp_trinity)")
    args = p.parse_args()
    args.skip = ["pytransrate" if s == "transrate" else s for s in args.skip]
    return args


def main():
    # Before anything prints: block-buffered stdout puts the banner and the
    # [step] lines out of order with the tools' output in a redirected log.
    line_buffer_stdio()
    Eval(parse_args()).main()


if __name__ == "__main__":
    main()
