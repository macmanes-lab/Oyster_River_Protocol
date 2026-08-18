#!/usr/bin/env python3
"""Python port of oyster.mk - the Oyster River Protocol assembly pipeline.

Usage:
    oyster.py --read1 R1.fq.gz --read2 R2.fq.gz --mem 110 --cpu 24 \\
              --runout myrun --strand RF

Runs the full pipeline end to end. Steps whose output files already exist
and are newer than their inputs are skipped, so a failed/interrupted run
can be re-invoked to resume where it left off.
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
import concurrent.futures
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
RED = "\033[31m"
RESET = "\033[0m"

SHORT_ASSEMBLY_NAMES = ("spades55.fasta", "spades75.fasta", "transabyss.fasta", "trinity.Trinity.fasta")

# Reference profile (minutes, from a representative run at --max-parallel 2)
# used only to decide submission order within the two remaining
# run_parallel() concurrent groups (orthofuser_branch vs. merge_branch;
# transrate vs. strandeval) -- run the historically slow step first so it
# isn't left waiting behind a quick one. The assemblers no longer go
# through run_parallel (see TRINITY_LANE_SHARE and the assembly lanes in
# main()), so they have no entries here. Each dataset is normally assembled
# only once, so this is a fixed relative ranking rather than something
# learned per-run. Names with no entry sort after every hinted step, in the
# order they were given.
STEP_TIME_HINTS = {
    "merge_branch": 27,
    "orthofuser_branch": 6,
    "transrate": 16,
    "strandeval": 2,
}

# Trinity's own wall time is overwhelmingly Phase 2 (thousands of small,
# independent per-gene-component assembly jobs dispatched via ParaFly --
# see run_trinity_phase2), which can only start once Phase 1 (Inchworm +
# Chrysalis prep, see run_trinity_phase1) has finished building the
# whole-transcriptome graph and partitioning reads. TRINITY_LANE_SHARE only
# governs that comparatively short Phase 1 window, split evenly with the
# short-assembler lane below -- Phase 1 is largely insensitive to its CPU
# share (Inchworm is hard-capped at --inchworm_cpu regardless of --CPU, and
# Chrysalis's clustering step is brief next to Phase 2). Phase 2 itself
# always claims the full --cpu/--mem budget once the short lane is done,
# rather than being capped at a fixed share for the rest of the run (see
# main()) -- real-run timing showed the short lane finishing in ~15h while
# Trinity's Phase 2 was only ~12% complete, i.e. days of the 20% share
# sitting idle under a fixed-for-the-whole-run split.
TRINITY_LANE_SHARE = 0.5

# A step that fails on a cluster is often transient (node preemption,
# filesystem hiccup, scheduler blip) rather than a real bug, so retry before
# giving up -- previously any single failed subprocess call killed the whole
# run immediately, discarding hours of unrelated concurrent work in the
# other lane. This does not help steps that fail deterministically (e.g. a
# missing dependency), which will just fail the same way on every attempt.
STEP_RETRIES = 2
STEP_RETRY_DELAY = 60


def awk_first_field(src: Path, dst: Path) -> None:
    with open(src) as inf, open(dst, "w") as outf:
        for line in inf:
            fields = line.split()
            outf.write((fields[0] if fields else "") + "\n")


def extract_gene_ids(diamond_txt: Path) -> set:
    ids = set()
    with open(diamond_txt) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            parts = cols[1].split("|")
            if len(parts) < 3:
                continue
            ids.add(parts[2].split("_")[0])
    return ids


def parse_unique_count(diamond_txt: Path) -> int:
    return len(extract_gene_ids(diamond_txt))


def write_sorted(path: Path, ids) -> None:
    with open(path, "w") as f:
        for i in sorted(ids):
            f.write(i + "\n")


def is_gzip(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def average_read_length(path: Path, n_records: int = 100) -> int:
    opener = gzip.open if is_gzip(path) else open
    lengths = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i >= n_records * 4:
                break
            if i % 4 == 1:
                lengths.append(len(line.strip()))
    return sum(lengths) // len(lengths) if lengths else 0


def hostname_suffix() -> str:
    parts = socket.gethostname().split(".")
    return ".".join(parts[2:5])


class Pipeline:
    def __init__(self, args):
        self.dir = Path(args.dir).resolve() if args.dir else Path.cwd()
        self.makedir = HERE
        self.cpu = args.cpu
        self.busco_threads = args.busco_threads or self.cpu
        self.mem = args.mem
        self.spades1_kmer = args.spades1_kmer
        self.spades2_kmer = args.spades2_kmer
        self.transabyss_kmer = args.transabyss_kmer
        self.read1 = Path(args.read1)
        self.read2 = Path(args.read2)
        self.runout = args.runout
        self.lineage = args.lineage
        self.strand = args.strand
        self.normalize_reads = args.normalize_reads
        self.tpm_filt = args.tpm_filt
        self.max_parallel = max(1, args.max_parallel)

        self.version = (self.makedir / "version.txt").read_text().strip()
        self.busco_config = self.makedir / "software" / "config.ini"
        os.environ["BUSCO_CONFIG_FILE"] = str(self.busco_config)
        self.diamond_db = self.makedir / "software" / "diamond" / "swissprot"

        self.reads_dir = self.dir / "reads"
        self.rcorr_dir = self.dir / "rcorr"
        self.assemblies_dir = self.dir / "assemblies"
        self.assemblies_working = self.assemblies_dir / "working"
        self.diamond_dir = self.assemblies_dir / "diamond"
        self.reports_dir = self.dir / "reports"
        self.orthofuse_dir = self.dir / "orthofuse" / self.runout
        self.orthofuse_working = self.orthofuse_dir / "working"
        self.quants_dir = self.dir / "quants"

        self.timing_log = self.reports_dir / f"{self.runout}.timing.log"
        self.run_cmd = "oyster.py " + " ".join(sys.argv[1:])
        self.steps = []
        self._timing_lock = threading.Lock()

    # -- process helpers -------------------------------------------------

    def run(self, cmd, cwd=None, retries=STEP_RETRIES, retry_delay=STEP_RETRY_DELAY,
            retry_cleanup=None, **kwargs):
        """retry_cleanup: path or iterable of paths to rmtree before each retry --
        for tools (SPAdes, TransAByss) that refuse to reuse a non-empty output
        dir rather than resuming, so a bare retry would just fail differently
        instead of actually re-attempting the work. Not needed for tools like
        Trinity that resume from their own checkpoints in-place.
        """
        printable = " ".join(str(c) for c in cmd)
        for attempt in range(retries + 1):
            print("+", printable)
            try:
                subprocess.run(cmd, check=True, cwd=str(cwd or self.dir), **kwargs)
                return
            except subprocess.CalledProcessError as e:
                if attempt == retries:
                    raise
                print(
                    f"*** step failed (exit {e.returncode}), retrying in {retry_delay}s "
                    f"[attempt {attempt + 2}/{retries + 1}] ***"
                )
                if retry_cleanup is not None:
                    paths = [retry_cleanup] if isinstance(retry_cleanup, (str, Path)) else retry_cleanup
                    for path in paths:
                        shutil.rmtree(path, ignore_errors=True)
                time.sleep(retry_delay)

    def conda_run(self, env, *cmd, **kwargs):
        self.run(["conda", "run", "--no-capture-output", "-n", env, *[str(c) for c in cmd]], **kwargs)

    # -- resumability ------------------------------------------------------

    def needs_run(self, outputs, inputs) -> bool:
        if not outputs:
            return True
        outputs = [Path(o) for o in outputs]
        if not all(o.exists() for o in outputs):
            return True
        existing_inputs = [Path(i) for i in inputs if i and Path(i).exists()]
        if not existing_inputs:
            return False
        newest_input = max(i.stat().st_mtime for i in existing_inputs)
        oldest_output = min(o.stat().st_mtime for o in outputs)
        return newest_input > oldest_output

    def step(self, name, outputs, inputs, func, timed=True):
        if not self.needs_run(outputs, inputs):
            print(f"[{name}] up to date, skipping")
            return
        print(f"\n=== {name} ===")
        start = time.time()
        func()
        elapsed = int(time.time() - start)
        if timed:
            self._record_timing(name, elapsed)

    def _record_timing(self, name, elapsed):
        with self._timing_lock:
            self.steps.append((name, elapsed))
            with open(self.timing_log, "a") as f:
                f.write(f"{name}\t{elapsed}\n")

    def run_parallel(self, jobs, max_workers=2):
        """Run independent (name, outputs, inputs, func) jobs concurrently.

        Each pending job's func is called as func(cpu=<split>, mem=<split>),
        where self.cpu/self.mem are divided across however many jobs actually
        run at once (capped at max_workers). Already up-to-date jobs are
        skipped and don't count toward the split.
        """
        pending = []
        for name, outputs, inputs, func in jobs:
            if self.needs_run(outputs, inputs):
                pending.append((name, func))
            else:
                print(f"[{name}] up to date, skipping")
        if not pending:
            return

        # Longest-processing-time-first: submit the slowest known jobs first so
        # they start immediately instead of waiting behind quick ones, using
        # a fixed reference profile (STEP_TIME_HINTS) rather than this run's
        # own history -- each dataset is normally only ever assembled once,
        # so there's no prior run of *this* data to learn from. Jobs with no
        # hint keep their given relative order, after every hinted job.
        pending.sort(key=lambda item: STEP_TIME_HINTS.get(item[0], -1), reverse=True)

        workers = min(max_workers, len(pending), max(1, self.cpu))
        job_cpu = max(1, self.cpu // workers)
        job_mem = max(1, self.mem // workers)
        if workers > 1:
            print(f"\n=== running {len(pending)} step(s), {workers} at a time (cpu={job_cpu}, mem={job_mem}G each) ===")

        def run_one(name, func):
            print(f"\n=== {name} (cpu={job_cpu}, mem={job_mem}G) ===")
            start = time.time()
            func(cpu=job_cpu, mem=job_mem)
            self._record_timing(name, int(time.time() - start))

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            futures = [ex.submit(run_one, name, func) for name, func in pending]
            for future in concurrent.futures.as_completed(futures):
                future.result()

    # -- setup / preflight -------------------------------------------------

    def setup(self):
        for d in (
            self.reads_dir, self.assemblies_dir, self.rcorr_dir, self.reports_dir,
            self.orthofuse_dir, self.quants_dir, self.diamond_dir, self.assemblies_working,
        ):
            d.mkdir(parents=True, exist_ok=True)

    def timing_init(self):
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        self.timing_log.write_text(f"Command: {self.run_cmd}\n\n")

    def which_in_env(self, env, binary):
        try:
            result = subprocess.run(
                ["conda", "run", "-n", env, "which", binary],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
        except FileNotFoundError:
            return None
        path = result.stdout.strip()
        if result.returncode == 0 and path and os.path.basename(path) == binary:
            return path
        return None

    def check(self):
        if self.which_in_env("orp", "salmon"):
            print("SALMON installed")
        else:
            sys.exit("*** SALMON is not installed, must fix ***")

        transrate_bin = self.makedir / "software" / "orp-transrate" / "transrate"
        if os.access(transrate_bin, os.X_OK):
            print("TRANSRATE installed")
        else:
            sys.exit("*** TRANSRATE is not installed, must fix ***")

        for env, binary, label in (
            ("orp", "seqtk", "SEQTK"),
            ("orp_busco", "busco", "BUSCO"),
            ("orp", "mcl", "MCL"),
            ("orp_spades", "rnaspades.py", "SPADES"),
            ("orp_trinity", "Trinity", "TRINITY"),
            ("orp", "trimmomatic", "TRIMMOMATIC"),
            ("orp_transabyss", "transabyss", "TRANSABYSS"),
            ("orp", "run_rcorrector.pl", "RCORRECTOR"),
            ("orp_orthofinder", "orthofinder", "ORTHOFINDER"),
        ):
            if self.which_in_env(env, binary):
                print(f"{label} installed")
            else:
                sys.exit(f"*** {label} is not installed, must fix ***")

    def welcome(self):
        print(RED)
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("                    _.-~~~~~-._")
        print("                 .-~   o   o   ~-.")
        print("                (   .-'~~~~~'-.   )        OYSTER RIVER PROTOCOL")
        print(f"                 '-.___________.-'         version {self.version}")
        print("                     '-.___.-'")
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + RESET + "\n")

    def readcheck(self):
        if not self.read1.exists():
            sys.exit("\n\n\n\n ERROR: YOUR READ1 FILE DOES NOT EXIST AT THE LOCATION YOU SPECIFIED\n\n\n\n ")
        if not self.read2.exists():
            sys.exit("\n\n\n\n ERROR: YOUR READ2 FILE DOES NOT EXIST AT THE LOCATION YOU SPECIFIED\n\n\n\n ")
        len1 = average_read_length(self.read1)
        len2 = average_read_length(self.read2)
        if not (len1 > self.spades2_kmer and len2 > self.spades2_kmer):
            sys.exit(
                f"\n\n\n\n IT LOOKS LIKE YOUR READS ARE NOT AT LEAST {self.spades2_kmer} BP LONG,\n "
                'PLEASE EDIT YOUR COMMAND USING THE "SPADES2_KMER=INT" FLAGS,\n'
                " SETTING THE ASSEMBLY KMER LENGTH TO AN ODD NUMBER LESS THAN YOUR READ LENGTH \n\n\n\n"
            )

    # -- trimming / correction ----------------------------------------------

    def trim1(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_1P.fastq"

    def trim2(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_2P.fastq"

    def cor1(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_1P.cor.fq"

    def cor2(self):
        return self.rcorr_dir / f"{self.runout}.TRIM_2P.cor.fq"

    def run_trimmomatic(self):
        baseout = self.rcorr_dir / f"{self.runout}.TRIM.fastq"
        clip = self.makedir / "barcodes" / "barcodes.fa"
        pe_args = [
            "PE", "-threads", str(self.cpu), "-baseout", str(baseout),
            str(self.read1), str(self.read2),
            "LEADING:3", "TRAILING:3",
            f"ILLUMINACLIP:{clip}:2:30:10:8:TRUE", "MINLEN:25",
        ]
        if hostname_suffix() == "bridges.psc.edu":
            trimmomatic_home = os.environ.get("TRIMMOMATIC_HOME", "")
            jar = f"{trimmomatic_home}/trimmomatic-0.36.jar"
            self.run(["java", f"-Xmx{self.mem}G", "-jar", jar, *pe_args])
        else:
            env = dict(os.environ, _JAVA_OPTIONS=f"-Xmx{self.mem}G")
            self.run(
                ["conda", "run", "--no-capture-output", "-n", "orp", "trimmomatic", *pe_args],
                env=env,
            )

    def run_rcorrector(self):
        self.conda_run(
            "orp", "run_rcorrector.pl", "-t", self.cpu, "-k", "31",
            "-1", self.trim1(), "-2", self.trim2(), "-od", self.rcorr_dir,
        )

    # -- assemblers ----------------------------------------------------------

    def trinity_out_dir(self):
        return self.assemblies_dir / f"{self.runout}.trinity"

    def _trinity_base_cmd(self, cpu, mem):
        cmd = ["Trinity"]
        if self.strand == "RF":
            cmd += ["--SS_lib_type", "RF"]
        elif self.strand == "FR":
            cmd += ["--SS_lib_type", "FR"]
        cmd += ["--no_version_check", "--bypass_java_version_check"]
        if not self.normalize_reads:
            cmd.append("--no_normalize_reads")
        cmd += [
            "--seqType", "fq",
            "--output", str(self.trinity_out_dir()),
            "--max_memory", f"{mem}G",
            "--left", str(self.cor1()), "--right", str(self.cor2()),
            "--CPU", str(cpu), "--inchworm_cpu", "10",
        ]
        return cmd

    def run_trinity_phase1(self, cpu=None, mem=None):
        # Inchworm + Chrysalis prep only, via Trinity's documented staged-
        # execution flag (https://github.com/trinityrnaseq/trinityrnaseq/wiki/
        # Running-Trinity#running-trinity-in-multiple-sequential-stages):
        # builds the whole-transcriptome contig graph, partitions reads per
        # gene component, writes the Phase-2 command list, then stops --
        # doesn't touch the actual per-component assembly (run_trinity_phase2
        # below). Runs alongside the short-assembler lane at TRINITY_LANE_SHARE
        # of --cpu/--mem.
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        cmd = self._trinity_base_cmd(cpu, mem) + ["--no_distributed_trinity_exec"]
        self.conda_run("orp_trinity", *cmd, retries=0)

    def run_trinity_phase2(self, cpu=None, mem=None):
        # Same command, no stop flag: per the docs above, Trinity resumes
        # from Phase 1's on-disk checkpoints straight into Phase 2 -- the
        # thousands of small, independent, single-threaded per-component
        # assembly jobs dispatched via ParaFly, and by far Trinity's
        # dominant cost. Runs alone (short-assembler lane has already
        # finished, see main()), so it gets the full machine rather than
        # splitting further.
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        out = self.assemblies_dir / f"{self.runout}.trinity.Trinity.fasta"
        cmd = self._trinity_base_cmd(cpu, mem) + ["--full_cleanup"]
        # No retries: this step's wall time dwarfs every other (hours to
        # days), so blindly retrying a deterministic failure could multiply
        # the wall time before finally giving up. It also resumes from its
        # own checkpoints in-place (same as Phase 1 above), so a manual
        # re-run of oyster.py after a transient failure loses little anyway.
        self.conda_run("orp_trinity", *cmd, retries=0)
        tmp = out.with_suffix(".fa")
        awk_first_field(out, tmp)
        tmp.replace(out)
        for f in self.assemblies_dir.glob("*gene_trans_map"):
            f.unlink()

    def run_spades(self, kmer, outname, workdir_suffix, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        mem = self.mem if mem is None else mem
        out = self.assemblies_dir / f"{self.runout}.{outname}.fasta"
        workdir = self.assemblies_dir / f"{self.runout}.spades_k{workdir_suffix}"
        cmd = ["rnaspades.py"]
        if self.strand == "RF":
            cmd.append("--ss-rf")
        elif self.strand == "FR":
            cmd.append("--ss-fr")
        cmd += [
            "--only-assembler", "-o", str(workdir),
            "--threads", str(cpu), "--memory", str(mem), "-k", str(kmer),
            "-1", str(self.cor1()), "-2", str(self.cor2()),
        ]
        self.conda_run("orp_spades", *cmd, retry_cleanup=workdir)
        shutil.move(str(workdir / "transcripts.fasta"), str(out))
        shutil.rmtree(workdir, ignore_errors=True)

    def run_spades55(self, cpu=None, mem=None):
        self.run_spades(self.spades1_kmer, "spades55", "55", cpu=cpu, mem=mem)

    def run_spades75(self, cpu=None, mem=None):
        self.run_spades(self.spades2_kmer, "spades75", "75", cpu=cpu, mem=mem)

    def run_transabyss(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        out = self.assemblies_dir / f"{self.runout}.transabyss.fasta"
        workdir = self.assemblies_dir / f"{self.runout}.transabyss"
        cmd = ["transabyss"]
        if self.strand in ("RF", "FR"):
            cmd.append("--SS")
        cmd += [
            "--threads", str(cpu), "--outdir", str(workdir),
            "--kmer", str(self.transabyss_kmer), "--length", "250",
            "--name", f"{self.runout}.transabyss.fasta",
            "--pe", str(self.cor1()), str(self.cor2()),
        ]
        self.conda_run("orp_transabyss", *cmd, retry_cleanup=workdir)
        final = workdir / f"{self.runout}.transabyss.fasta-final.fa"
        awk_first_field(final, out)
        shutil.rmtree(workdir, ignore_errors=True)

    # -- orthofuse merge -------------------------------------------------------

    def short_fasta_paths(self):
        return [self.orthofuse_working / f"{self.runout}.{n}.short.fasta" for n in SHORT_ASSEMBLY_NAMES]

    def run_filtershort(self):
        self.orthofuse_working.mkdir(parents=True, exist_ok=True)
        for n in ("transabyss.fasta", "spades75.fasta", "spades55.fasta", "trinity.Trinity.fasta"):
            fasta = self.assemblies_dir / f"{self.runout}.{n}"
            outp = self.orthofuse_working / f"{fasta.name}.short.fasta"
            self.conda_run("orp", "python", self.makedir / "scripts" / "long.seq.py", fasta, outp, "200")

    def run_orthofuser(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        self.conda_run(
            "orp_orthofinder", "orthofinder",
            "-d", "-I", "12", "-f", self.orthofuse_working, "-og", "-t", cpu, "-a", cpu,
        )
        (self.orthofuse_dir / "orthofuser.done").touch()

    def merge(self):
        out = self.orthofuse_dir / "merged.fasta"
        with open(out, "wb") as outf:
            for p in self.short_fasta_paths():
                with open(p, "rb") as inf:
                    shutil.copyfileobj(inf, outf)

    def find_orthogroups_txt(self):
        match = next(self.orthofuse_working.rglob("Orthogroups.txt"), None)
        if match is None:
            sys.exit("Orthogroups.txt not found under orthofuse working directory")
        return match

    def makelist(self):
        orthogroups = self.find_orthogroups_txt()
        with open(orthogroups) as f:
            n = sum(1 for _ in f)
        out = self.orthofuse_dir / f"{self.runout}.list"
        with open(out, "w") as f:
            for i in range(1, n + 1):
                f.write(f"{i}\n")

    def makegroups(self):
        orthogroups = self.find_orthogroups_txt()
        with open(orthogroups) as f:
            lines = f.readlines()

        def write_group(i):
            tokens = lines[i - 1].split()
            group_file = self.orthofuse_dir / f"{i}.groups"
            with open(group_file, "w") as f:
                for tok in tokens[1:]:
                    f.write(tok + "\n")

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpu) as ex:
            list(ex.map(write_group, range(1, len(lines) + 1)))
        (self.orthofuse_dir / "groups.done").touch()

    def orthotransrate(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        outdir = self.orthofuse_dir / "merged"
        self.conda_run(
            "orp", self.makedir / "software" / "orp-transrate" / "transrate",
            "-o", outdir, "-t", cpu, "-a", self.orthofuse_dir / "merged.fasta",
            "--left", self.cor1(), "--right", self.cor2(),
        )
        for f in outdir.rglob("*.bam"):
            f.unlink()

    def makeorthout(self):
        print("Picking the best contig per orthogroup")
        contigs_csv = next(self.orthofuse_dir.rglob("contigs.csv"), None)
        if contigs_csv is None:
            sys.exit("contigs.csv not found under orthofuse directory")
        good_list = self.orthofuse_dir / f"good.{self.runout}.list"
        self.conda_run(
            "orp", "python", self.makedir / "scripts" / "pick_best_contigs.py",
            contigs_csv, self.orthofuse_dir, good_list,
        )
        for f in self.orthofuse_dir.glob("*groups"):
            f.unlink()

    def orthofusing(self):
        good_list = self.orthofuse_dir / f"good.{self.runout}.list"
        out = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        with open(out, "w") as outf:
            subprocess.run(
                ["conda", "run", "--no-capture-output", "-n", "orp", "python",
                 str(self.makedir / "scripts" / "filter.py"), str(self.orthofuse_dir / "merged.fasta"), str(good_list)],
                check=True, stdout=outf, cwd=self.dir,
            )

    # -- diamond passes ----------------------------------------------------

    def diamond_jobs(self):
        return [
            (self.assemblies_dir / f"{self.runout}.orthomerged.fasta", self.diamond_dir / f"{self.runout}.orthomerged.diamond.txt"),
            (self.assemblies_dir / f"{self.runout}.transabyss.fasta", self.diamond_dir / f"{self.runout}.transabyss.diamond.txt"),
            (self.assemblies_dir / f"{self.runout}.spades75.fasta", self.diamond_dir / f"{self.runout}.spades75.diamond.txt"),
            (self.assemblies_dir / f"{self.runout}.spades55.fasta", self.diamond_dir / f"{self.runout}.spades55.diamond.txt"),
            (self.assemblies_dir / f"{self.runout}.trinity.Trinity.fasta", self.diamond_dir / f"{self.runout}.trinity.diamond.txt"),
        ]

    def run_diamond_one(self, query, out, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        self.conda_run(
            "orp", "diamond", "blastx", "--quiet", "-p", cpu,
            "-e", "1e-8", "--top", "0.1", "-q", query, "-d", self.diamond_db, "-o", out,
        )

    def diamond_uniq(self):
        mapping = {
            "trinity": self.diamond_dir / f"{self.runout}.trinity.diamond.txt",
            "sp75": self.diamond_dir / f"{self.runout}.spades75.diamond.txt",
            "sp55": self.diamond_dir / f"{self.runout}.spades55.diamond.txt",
            "transabyss": self.diamond_dir / f"{self.runout}.transabyss.diamond.txt",
        }
        for label, path in mapping.items():
            count = parse_unique_count(path)
            (self.diamond_dir / f"{self.runout}.unique.{label}.txt").write_text(f"{count}\n")

    def make_list1(self):
        ids = extract_gene_ids(self.diamond_dir / f"{self.runout}.orthomerged.diamond.txt")
        write_sorted(self.diamond_dir / f"{self.runout}.list1", ids)

    def make_list2(self):
        ids = set()
        for name in ("transabyss", "spades75", "spades55", "trinity"):
            ids |= extract_gene_ids(self.diamond_dir / f"{self.runout}.{name}.diamond.txt")
        write_sorted(self.diamond_dir / f"{self.runout}.list2", ids)

    def make_list3(self):
        list1 = set(line.rstrip("\n") for line in open(self.diamond_dir / f"{self.runout}.list1"))
        list2_path = self.diamond_dir / f"{self.runout}.list2"
        out = self.diamond_dir / f"{self.runout}.list3"
        with open(list2_path) as f, open(out, "w") as o:
            for line in f:
                if line.rstrip("\n") not in list1:
                    o.write(line)

    def make_list5(self):
        self.conda_run(
            "orp", "python", self.makedir / "scripts" / "build_list5.py",
            self.diamond_dir / f"{self.runout}.list3", self.diamond_dir / f"{self.runout}.list5",
            self.diamond_dir / f"{self.runout}.transabyss.diamond.txt",
            self.diamond_dir / f"{self.runout}.spades75.diamond.txt",
            self.diamond_dir / f"{self.runout}.spades55.diamond.txt",
            self.diamond_dir / f"{self.runout}.trinity.diamond.txt",
        )

    def make_list6(self):
        fasta = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        out = self.diamond_dir / f"{self.runout}.list6"
        with open(fasta) as f, open(out, "w") as o:
            for line in f:
                if line.startswith(">"):
                    o.write(line[1:])

    def make_list7(self):
        list6 = set(line.rstrip("\n") for line in open(self.diamond_dir / f"{self.runout}.list6"))
        list5_path = self.diamond_dir / f"{self.runout}.list5"
        out = self.diamond_dir / f"{self.runout}.list7"
        with open(list5_path) as f, open(out, "w") as o:
            for line in f:
                if line.rstrip("\n") not in list6:
                    o.write(line)

    def posthack(self):
        fastas = " ".join(
            str(self.assemblies_dir / f"{self.runout}.{n}")
            for n in ("spades55.fasta", "spades75.fasta", "transabyss.fasta", "trinity.Trinity.fasta")
        )
        list7 = self.diamond_dir / f"{self.runout}.list7"
        newbies = self.diamond_dir / f"{self.runout}.newbies.fasta"
        orthomerged = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        working_out = self.assemblies_working / f"{self.runout}.orthomerged.fasta"
        filter_py = self.makedir / "scripts" / "filter.py"
        script = f"python {filter_py} <(cat {fastas}) {list7} >> {newbies}"
        self.run(["conda", "run", "--no-capture-output", "-n", "orp", "bash", "-c", script])
        with open(working_out, "wb") as outf:
            for p in (newbies, orthomerged):
                with open(p, "rb") as inf:
                    shutil.copyfileobj(inf, outf)

    # -- dedup / quantify -----------------------------------------------------

    def cdhit(self):
        src = self.assemblies_working / f"{self.runout}.orthomerged.fasta"
        out = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        self.conda_run(
            "orp", "cd-hit-est", "-M", self.mem * 1000, "-T", self.cpu,
            "-c", ".98", "-i", src, "-o", out,
        )

    def orp_diamond(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        src = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        out = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        self.conda_run(
            "orp", "diamond", "blastx", "--quiet", "-p", cpu,
            "-e", "1e-8", "--top", "0.1", "-q", src, "-d", self.diamond_db, "-o", out,
        )
        out.touch()

    def orp_uniq(self):
        diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        count = parse_unique_count(diamond_txt)
        (self.assemblies_working / f"{self.runout}.unique.ORP.txt").write_text(f"{count}\n")
        (self.assemblies_working / f"{self.runout}.unique.ORP.done").touch()

    def salmon_index(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        src = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        idx = self.quants_dir / f"{self.runout}.ortho.idx"
        self.conda_run(
            "orp", "salmon", "index", "--no-version-check", "-t", src,
            "-i", idx, "-k", "31", "--threads", cpu,
        )

    def salmon(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        idx = self.quants_dir / f"{self.runout}.ortho.idx"
        outdir = self.quants_dir / f"salmon_orthomerged_{self.runout}"
        self.conda_run(
            "orp", "salmon", "quant", "--no-version-check", "--validateMappings",
            "-p", cpu, "-i", idx, "--seqBias", "--gcBias", "--libType", "A",
            "-1", self.cor1(), "-2", self.cor2(), "-o", outdir,
        )

    def filter_tpm(self):
        quant = self.quants_dir / f"salmon_orthomerged_{self.runout}" / "quant.sf"
        high = self.assemblies_working / f"{self.runout}.HIGHEXP.txt"
        low = self.assemblies_working / f"{self.runout}.LOWEXP.txt"
        with open(quant) as f, open(high, "w") as hf, open(low, "w") as lf:
            next(f, None)
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4:
                    continue
                try:
                    tpm = float(cols[3])
                except ValueError:
                    continue
                if tpm > self.tpm_filt:
                    hf.write(cols[0] + "\n")
                elif tpm < self.tpm_filt:
                    lf.write(cols[0] + "\n")
        (self.assemblies_dir / f"{self.runout}.filter.done").touch()
        print("\n\n\n\n PART: TPM_FILT MAKE LOW AND HIGH\n\n\n\n")

    def secondfilter(self):
        intermediate = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        low = self.assemblies_working / f"{self.runout}.LOWEXP.txt"
        high = self.assemblies_working / f"{self.runout}.HIGHEXP.txt"
        diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        filter_py = self.makedir / "scripts" / "filter.py"

        if low.exists() and low.stat().st_size > 0:
            before = self.assemblies_working / f"{self.runout}.ORP_BEFORE_TPM_FILT.fasta"
            shutil.copy(intermediate, before)

            highexp_fasta = self.assemblies_working / f"{self.runout}.ORP.HIGHEXP.fasta"
            with open(highexp_fasta, "w") as outf:
                subprocess.run(
                    ["conda", "run", "--no-capture-output", "-n", "orp", "python",
                     str(filter_py), str(intermediate), str(high)],
                    check=True, stdout=outf, cwd=self.dir,
                )

            low_ids = set(line.strip() for line in open(low) if line.strip())
            blasted = self.assemblies_working / f"{self.runout}.blasted"
            donotremove = self.assemblies_working / f"{self.runout}.donotremove.list"
            do_not_remove_ids = set()
            with open(diamond_txt) as f, open(blasted, "a") as bf:
                for line in f:
                    cols = line.rstrip("\n").split("\t")
                    if cols and cols[0] in low_ids:
                        bf.write(line)
                        do_not_remove_ids.add(cols[0])
            with open(donotremove, "a") as df:
                for i in sorted(do_not_remove_ids):
                    print(i)
                    df.write(i + "\n")

            saveme = self.assemblies_working / f"{self.runout}.saveme.fasta"
            with open(saveme, "w") as outf:
                subprocess.run(
                    ["conda", "run", "--no-capture-output", "-n", "orp", "python",
                     str(filter_py), str(intermediate), str(donotremove)],
                    check=True, stdout=outf, cwd=self.dir,
                )

            with open(orp_fasta, "wb") as outf:
                for p in (saveme, highexp_fasta):
                    with open(p, "rb") as inf:
                        shutil.copyfileobj(inf, outf)
            print("\n\n\n\n PART: FILTER LOW STUFF \n\n\n\n")
        else:
            shutil.copy(intermediate, orp_fasta)
            print("\n\n\n\n PART: THERE IS NO LOW STUFF \n\n\n\n")

    # -- QC / reporting ----------------------------------------------------

    def busco(self, cpu=None, mem=None):
        cpu = self.busco_threads if cpu is None else cpu
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        self.conda_run(
            "orp_busco", "busco", "--offline", "--lineage", self.lineage,
            "--download_path", self.makedir / "busco_dbs",
            "-i", orp_fasta, "-m", "transcriptome", "--cpu", cpu,
            "-o", f"run_{self.runout}.ORP", "--config", self.busco_config,
        )
        for p in self.dir.glob(f"run_{self.runout}*"):
            shutil.move(str(p), str(self.reports_dir / p.name))
        (self.reports_dir / f"{self.runout}.busco.done").touch()

    def transrate(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        outdir = self.reports_dir / f"transrate_{self.runout}"
        self.conda_run(
            "orp", self.makedir / "software" / "orp-transrate" / "transrate",
            "-o", outdir, "-a", orp_fasta,
            "--left", self.cor1(), "--right", self.cor2(), "-t", cpu,
        )
        for f in outdir.rglob("*.bam"):
            f.unlink()

    def trinity_perllib_dir(self):
        result = subprocess.run(
            ["conda", "run", "-n", "orp_trinity", "which", "Trinity"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True,
        )
        trinity_path = Path(result.stdout.strip()).resolve()
        return trinity_path.parent / "PerlLib"

    def strandeval(self, cpu=None, mem=None):
        cpu = self.cpu if cpu is None else cpu
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        r1, r2 = self.cor1(), self.cor2()
        self.conda_run("orp_trinity", "bwa", "index", "-p", self.runout, orp_fasta)

        pipeline_script = (
            f'conda run --no-capture-output -n orp_trinity bash -c '
            f'"bwa mem -t {cpu} {self.runout} '
            f'<(seqtk sample -s 23894 {r1} 400000) <(seqtk sample -s 23894 {r2} 400000)" '
            f"| conda run --no-capture-output -n orp samtools view -@10 -Sb - "
            f"| conda run --no-capture-output -n orp samtools sort -T {self.runout} -O bam -@10 "
            f"-o {self.runout}.sorted.bam -"
        )
        self.run(["bash", "-o", "pipefail", "-c", pipeline_script])

        flagstat_out = self.assemblies_dir / f"{self.runout}.flagstat"
        with open(flagstat_out, "w") as outf:
            subprocess.run(
                ["conda", "run", "--no-capture-output", "-n", "orp", "samtools", "flagstat",
                 f"{self.runout}.sorted.bam"],
                check=True, stdout=outf, cwd=self.dir,
            )

        perllib = self.trinity_perllib_dir()
        self.conda_run(
            "orp_trinity", "perl", "-I", str(perllib),
            str(self.makedir / "scripts" / "examine_strand.pl"),
            f"{self.runout}.sorted.bam", self.runout,
        )

        dat_file = self.dir / f"{self.runout}.dat"
        hist_input = self.dir / f"{self.runout}.hist_input.txt"
        with open(dat_file) as f, open(hist_input, "w") as out:
            next(f, None)
            for line in f:
                cols = line.rstrip("\n").split()
                if len(cols) >= 5:
                    out.write(cols[4] + "\n")

        hist_result = subprocess.run(
            ["conda", "run", "--no-capture-output", "-n", "orp_trinity", "bash", "-c",
             f"hist -p '#' -c red {hist_input}"],
            check=True, cwd=self.dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
        )

        if hist_input.exists():
            hist_input.unlink()
        sorted_bam = self.dir / f"{self.runout}.sorted.bam"
        if sorted_bam.exists():
            sorted_bam.unlink()
        for ext in ("bwt", "pac", "ann", "amb", "sa", "dat"):
            p = self.dir / f"{self.runout}.{ext}"
            if p.exists():
                p.unlink()

        (self.reports_dir / f"{self.runout}.strandeval.done").touch()
        histogram_text = hist_result.stdout.rstrip("\n")
        summary = (
            "\n*****  STRAND EXAMINATION HISTOGRAM ***** \n"
            f"{histogram_text}\n"
            "\n*****  See the following link for interpretation ***** \n"
            "*****  https://oyster-river-protocol.readthedocs.io/en/latest/strandexamine.html ***** \n"
        )
        print(summary)
        (self.reports_dir / f"{self.runout}.strandeval_summary.txt").write_text(summary)

    def reportgen(self):
        runout = self.runout
        lines = []

        def emit(label, value):
            text = f"{label}      {value}"
            print(text)
            lines.append(text)

        print(f"\n\n*****  QUALITY REPORT FOR: {runout} using the ORP version {self.version} ****")
        orp_fasta = self.assemblies_dir / f"{runout}.ORP.fasta"
        print(f"\n*****  THE ASSEMBLY CAN BE FOUND HERE: {orp_fasta} **** \n")

        busco_line = ""
        busco_dir = self.reports_dir / f"run_{runout}.ORP"
        for p in busco_dir.rglob("short*.txt"):
            for line in open(p):
                if re.match(r"^\s*C:[0-9]", line):
                    busco_line = line.strip()
        emit("*****  BUSCO SCORE ~~~~~~~~~~~~~~~~~~~~~~>", busco_line)

        csv_path = next((self.reports_dir / f"transrate_{runout}").rglob("assemblies.csv"), None)
        rows = list(csv.reader(open(csv_path))) if csv_path else []
        transrate_score = rows[1][36] if len(rows) > 1 and len(rows[1]) > 36 else ""
        transrate_optimal = rows[1][37] if len(rows) > 1 and len(rows[1]) > 37 else ""
        emit("*****  TRANSRATE SCORE ~~~~~~~~~~~~~~~~~~>     ", transrate_score)
        emit("*****  TRANSRATE OPTIMAL SCORE ~~~~~~~~~~>     ", transrate_optimal)

        def read_count(path):
            return path.read_text().strip() if path.exists() else ""

        emit("*****  UNIQUE GENES ORP ~~~~~~~~~~~~~~~~~>     ", read_count(self.assemblies_working / f"{runout}.unique.ORP.txt"))
        emit("*****  UNIQUE GENES TRINITY ~~~~~~~~~~~~~>     ", read_count(self.diamond_dir / f"{runout}.unique.trinity.txt"))
        emit("*****  UNIQUE GENES SPADES55 ~~~~~~~~~~~~>     ", read_count(self.diamond_dir / f"{runout}.unique.sp55.txt"))
        emit("*****  UNIQUE GENES SPADES75 ~~~~~~~~~~~~>     ", read_count(self.diamond_dir / f"{runout}.unique.sp75.txt"))
        emit("*****  UNIQUE GENES TRANSABYSS ~~~~~~~~~~>     ", read_count(self.diamond_dir / f"{runout}.unique.transabyss.txt"))

        proper_pairs = ""
        flagstat = self.assemblies_dir / f"{runout}.flagstat"
        if flagstat.exists():
            for line in open(flagstat):
                if "properly paired" in line:
                    fields = line.split()
                    if len(fields) > 5:
                        proper_pairs = fields[5].lstrip("(")
        emit("*****  READS MAPPED AS PROPER PAIRS ~~~~~>     ", proper_pairs)
        print(" \n")

        strandeval_summary_path = self.reports_dir / f"{runout}.strandeval_summary.txt"
        if strandeval_summary_path.exists():
            strandeval_summary = strandeval_summary_path.read_text()
            print(strandeval_summary)
            lines.append(strandeval_summary.rstrip("\n"))

        qualreport_path = self.reports_dir / f"qualreport.{runout}"
        qualreport_path.write_text("\n".join(lines) + "\n")
        (self.reports_dir / f"qualreport.{runout}.done").touch()

    def timing_report(self, wall_clock_seconds):
        with open(self.timing_log) as f:
            all_lines = f.readlines()
        header = all_lines[0].rstrip("\n") if all_lines else ""
        body = [l for l in all_lines if "\t" in l]

        out_lines = [header, "", f"*****  STEP TIMING for {self.runout} ***** ", ""]
        for line in body:
            name, secs = line.rstrip("\n").split("\t")
            secs = int(secs)
            h, m, s = secs // 3600, (secs % 3600) // 60, secs % 60
            out_lines.append(f"{name:<16} {h:02d}:{m:02d}:{s:02d}")
        # TOTAL is real wall-clock elapsed time, not a sum of the lines above:
        # concurrent jobs overlap (their durations shouldn't add), and a
        # branch's own entry already includes the sub-steps it calls, so
        # summing every line double-counts and overstates the true runtime.
        h, m, s = wall_clock_seconds // 3600, (wall_clock_seconds % 3600) // 60, wall_clock_seconds % 60
        out_lines.append(f"{'TOTAL':<16} {h:02d}:{m:02d}:{s:02d}")

        text = "\n".join(out_lines) + "\n"
        print(text)
        self.timing_log.write_text(text)
        print(f"\nFull timing log saved to: {self.timing_log}\n")

    # -- orchestration -------------------------------------------------------

    def main(self):
        pipeline_start = time.time()
        self.setup()
        self.timing_init()
        self.check()
        self.welcome()
        self.readcheck()

        t1, t2 = self.trim1(), self.trim2()
        c1, c2 = self.cor1(), self.cor2()
        trinity_fa = self.assemblies_dir / f"{self.runout}.trinity.Trinity.fasta"
        sp75 = self.assemblies_dir / f"{self.runout}.spades75.fasta"
        sp55 = self.assemblies_dir / f"{self.runout}.spades55.fasta"
        ta = self.assemblies_dir / f"{self.runout}.transabyss.fasta"
        short_fastas = self.short_fasta_paths()
        orthofuser_done = self.orthofuse_dir / "orthofuser.done"
        merged_fasta = self.orthofuse_dir / "merged.fasta"
        list_file = self.orthofuse_dir / f"{self.runout}.list"
        groups_done = self.orthofuse_dir / "groups.done"
        merged_csv = self.orthofuse_dir / "merged" / "assemblies.csv"
        good_list = self.orthofuse_dir / f"good.{self.runout}.list"
        orthomerged_fasta = self.assemblies_dir / f"{self.runout}.orthomerged.fasta"
        diamond_outs = [o for _, o in self.diamond_jobs()]
        diamond_orthomerged = self.diamond_dir / f"{self.runout}.orthomerged.diamond.txt"
        diamond_ta = self.diamond_dir / f"{self.runout}.transabyss.diamond.txt"
        diamond_sp75 = self.diamond_dir / f"{self.runout}.spades75.diamond.txt"
        diamond_sp55 = self.diamond_dir / f"{self.runout}.spades55.diamond.txt"
        diamond_trinity = self.diamond_dir / f"{self.runout}.trinity.diamond.txt"
        uniq_outs = [
            self.diamond_dir / f"{self.runout}.unique.trinity.txt",
            self.diamond_dir / f"{self.runout}.unique.sp75.txt",
            self.diamond_dir / f"{self.runout}.unique.sp55.txt",
            self.diamond_dir / f"{self.runout}.unique.transabyss.txt",
        ]
        list1 = self.diamond_dir / f"{self.runout}.list1"
        list2 = self.diamond_dir / f"{self.runout}.list2"
        list3 = self.diamond_dir / f"{self.runout}.list3"
        list5 = self.diamond_dir / f"{self.runout}.list5"
        list6 = self.diamond_dir / f"{self.runout}.list6"
        list7 = self.diamond_dir / f"{self.runout}.list7"
        newbies = self.diamond_dir / f"{self.runout}.newbies.fasta"
        working_orthomerged = self.assemblies_working / f"{self.runout}.orthomerged.fasta"
        orp_intermediate = self.assemblies_dir / f"{self.runout}.ORP.intermediate.fasta"
        orp_diamond_txt = self.assemblies_dir / f"{self.runout}.ORP.diamond.txt"
        unique_orp_done = self.assemblies_working / f"{self.runout}.unique.ORP.done"
        ortho_idx = self.quants_dir / f"{self.runout}.ortho.idx"
        quant_sf = self.quants_dir / f"salmon_orthomerged_{self.runout}" / "quant.sf"
        filter_done = self.assemblies_dir / f"{self.runout}.filter.done"
        low_txt = self.assemblies_working / f"{self.runout}.LOWEXP.txt"
        high_txt = self.assemblies_working / f"{self.runout}.HIGHEXP.txt"
        orp_fasta = self.assemblies_dir / f"{self.runout}.ORP.fasta"
        busco_done = self.reports_dir / f"{self.runout}.busco.done"
        transrate_csv = self.reports_dir / f"transrate_{self.runout}" / "assemblies.csv"
        strandeval_done = self.reports_dir / f"{self.runout}.strandeval.done"
        qualreport_done = self.reports_dir / f"qualreport.{self.runout}.done"

        self.step("run_trimmomatic", [t1, t2], [self.read1, self.read2], self.run_trimmomatic)
        self.step("run_rcorrector", [c1, c2], [t1, t2], self.run_rcorrector)

        # Two resource lanes rather than an even run_parallel() split across
        # all four assemblers -- see TRINITY_LANE_SHARE above for why.
        lane2_cpu = max(1, round(self.cpu * (1 - TRINITY_LANE_SHARE)))
        lane2_mem = max(1, round(self.mem * (1 - TRINITY_LANE_SHARE)))
        trinity_cpu = max(1, self.cpu - lane2_cpu)
        trinity_mem = max(1, self.mem - lane2_mem)

        def _lane_failed(lane_name, e):
            # ThreadPoolExecutor.__exit__ (below) calls shutdown(wait=True) even
            # while an exception unwinds, so the *other* lane's still-running
            # future keeps blocking the run's exit -- Trinity in particular can
            # still be tens of hours out. Log the failure the moment it happens
            # rather than leaving it silent until the other lane finally finishes.
            print(
                f"\n*** [{lane_name}] lane failed ({e}) -- other lane keeps running; "
                "this run will still abort once it finishes ***",
                flush=True,
            )

        def trinity_phase1_lane():
            try:
                self.step(
                    "run_trinity_phase1", [self.trinity_out_dir() / "recursive_trinity.cmds.ok"], [c1],
                    partial(self.run_trinity_phase1, cpu=trinity_cpu, mem=trinity_mem),
                )
            except Exception as e:
                _lane_failed("run_trinity_phase1", e)
                raise

        def short_assembler_lane():
            # Slowest of the three first; diamond_{transabyss,spades75,
            # spades55} depend only on their own assembly (not on Trinity or
            # the orthofuser merge below), so each fires as soon as its
            # assembly is done instead of waiting for the merge stage.
            try:
                for step_name, outputs, inputs, assemble, diamond_name, diamond_out in (
                    ("run_transabyss", [ta], [c1, c2], self.run_transabyss, "transabyss", diamond_ta),
                    ("run_spades55", [sp55], [c1, c2], self.run_spades55, "spades55", diamond_sp55),
                    ("run_spades75", [sp75], [c1, c2], self.run_spades75, "spades75", diamond_sp75),
                ):
                    self.step(step_name, outputs, inputs, partial(assemble, cpu=lane2_cpu, mem=lane2_mem))
                    query = outputs[0]
                    self.step(
                        f"diamond_{diamond_name}", [diamond_out], [query],
                        partial(self.run_diamond_one, query, diamond_out, cpu=lane2_cpu),
                    )
            except Exception as e:
                _lane_failed("short-assembler", e)
                raise

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
            for f in concurrent.futures.as_completed([ex.submit(trinity_phase1_lane), ex.submit(short_assembler_lane)]):
                f.result()

        # Both lanes have released the machine -- Trinity Phase 2 gets the
        # full --cpu/--mem budget, uncontested (see TRINITY_LANE_SHARE and
        # run_trinity_phase2 above).
        self.step(
            "run_trinity_phase2", [trinity_fa],
            [self.trinity_out_dir() / "recursive_trinity.cmds.ok"],
            self.run_trinity_phase2,
        )

        self.step("run_filtershort", short_fastas, [ta, sp75, sp55, trinity_fa], self.run_filtershort)

        def orthofuser_branch(cpu=None, mem=None):
            self.step("run_orthofuser", [orthofuser_done], short_fastas, partial(self.run_orthofuser, cpu=cpu))
            self.step("makelist", [list_file], [orthofuser_done], self.makelist, timed=False)
            self.step("makegroups", [groups_done], [list_file], self.makegroups, timed=False)

        def merge_branch(cpu=None, mem=None):
            self.step("merge", [merged_fasta], short_fastas, self.merge, timed=False)
            self.step("orthotransrate", [merged_csv], [merged_fasta, c1, c2], partial(self.orthotransrate, cpu=cpu))

        # run_orthofuser->makelist->makegroups and merge->orthotransrate are
        # independent chains that both only need short_fastas; they join at
        # makeorthout below.
        self.run_parallel(
            [
                ("orthofuser_branch", [groups_done], short_fastas, orthofuser_branch),
                ("merge_branch", [merged_csv], short_fastas, merge_branch),
            ],
            max_workers=self.max_parallel,
        )
        self.step("makeorthout", [good_list], [groups_done, merged_csv], self.makeorthout)
        self.step("orthofusing", [orthomerged_fasta], [good_list, merged_fasta], self.orthofusing)

        # diamond_transabyss/spades75/spades55 already ran in the short
        # assembler lane above. Only these two remain: orthomerged depends
        # on the merge stage just above, and trinity depends on the trinity
        # lane -- both are only just now guaranteed to be ready.
        print("\n\n\n\n Starting diamond \n\n\n\n")
        self.step(
            "diamond_orthomerged", [diamond_orthomerged], [orthomerged_fasta],
            partial(self.run_diamond_one, orthomerged_fasta, diamond_orthomerged),
        )
        self.step(
            "diamond_trinity", [diamond_trinity], [trinity_fa],
            partial(self.run_diamond_one, trinity_fa, diamond_trinity),
        )
        self.step("diamond_uniq", uniq_outs, diamond_outs, self.diamond_uniq, timed=False)
        self.step("make_list1", [list1], [diamond_orthomerged], self.make_list1, timed=False)
        self.step("make_list2", [list2], [diamond_trinity, diamond_sp75, diamond_sp55, diamond_ta], self.make_list2, timed=False)
        self.step("make_list3", [list3], [list1, list2], self.make_list3, timed=False)
        self.step("make_list5", [list5], [list3, diamond_ta, diamond_sp75, diamond_sp55, diamond_trinity], self.make_list5)
        self.step("make_list6", [list6], [orthomerged_fasta], self.make_list6, timed=False)
        self.step("make_list7", [list7], [list6, list5], self.make_list7, timed=False)
        self.step("posthack", [newbies, working_orthomerged], [list7], self.posthack)
        self.step("cdhit", [orp_intermediate], [working_orthomerged], self.cdhit)

        # orp_diamond is the same CPU-bound diamond blastx as above, paired
        # here with salmon_branch which is tiny (~2s); halving orp_diamond's
        # CPU to overlap with it costs more than the overlap saves, so both
        # run sequentially at full CPU instead.
        self.step("orp_diamond", [orp_diamond_txt], [orp_intermediate], self.orp_diamond)
        self.step("orp_uniq", [unique_orp_done], [orp_diamond_txt], self.orp_uniq, timed=False)
        self.step("salmon_index", [ortho_idx], [orp_intermediate], self.salmon_index)
        self.step("salmon", [quant_sf], [ortho_idx, c1, c2], self.salmon)
        self.step("filter", [filter_done], [orp_intermediate, quant_sf, orp_diamond_txt], self.filter_tpm, timed=False)
        self.step(
            "secondfilter", [orp_fasta],
            [filter_done, low_txt, high_txt, orp_intermediate, quant_sf, orp_diamond_txt],
            self.secondfilter,
        )
        # BUSCO dwarfs transrate/strandeval (minutes vs. seconds); splitting its
        # CPUs to overlap with them would slow it down far more than the
        # overlap could ever save, so it gets the full --cpu/--busco-threads
        # budget to itself. transrate and strandeval are independent of each
        # other and cheap, so they still run as a small concurrent pair.
        self.step("busco", [busco_done], [orp_fasta], self.busco)
        self.run_parallel(
            [
                ("transrate", [transrate_csv], [orp_fasta, c1, c2], self.transrate),
                ("strandeval", [strandeval_done], [orp_fasta], self.strandeval),
            ],
            max_workers=self.max_parallel,
        )
        self.step("reportgen", [qualreport_done], [unique_orp_done, orp_fasta], self.reportgen, timed=False)

        self.timing_report(int(time.time() - pipeline_start))


def parse_args():
    p = argparse.ArgumentParser(description="Python port of oyster.mk - the Oyster River Protocol pipeline.")
    version = (HERE / "version.txt").read_text().strip()
    p.add_argument("--version", action="version", version=f"Oyster River Protocol {version}")
    p.add_argument("--read1", required=True, help="path to R1 fastq(.gz)")
    p.add_argument("--read2", required=True, help="path to R2 fastq(.gz)")
    p.add_argument("--mem", type=int, default=110, help="memory in GB (default: 110)")
    p.add_argument("--cpu", type=int, default=16, help="CPU threads (default: 16)")
    p.add_argument("--busco-threads", type=int, default=None, help="BUSCO threads (default: same as --cpu)")
    p.add_argument("--runout", default="USER_RUN", help="run name prefix (default: USER_RUN)")
    p.add_argument("--strand", choices=["RF", "FR", ""], default="", help="strand-specificity (default: unstranded)")
    p.add_argument("--lineage", default="eukaryota_odb12.2", help="BUSCO lineage (default: eukaryota_odb12.2)")
    p.add_argument("--normalize-reads", action="store_true", help="let Trinity normalize reads (default: off, i.e. --no_normalize_reads)")
    p.add_argument("--tpm-filt", type=float, default=0, help="TPM filter threshold (default: 0)")
    p.add_argument("--spades1-kmer", type=int, default=55, help="rnaSPAdes k-mer for spades55 (default: 55)")
    p.add_argument("--spades2-kmer", type=int, default=75, help="rnaSPAdes k-mer for spades75 (default: 75)")
    p.add_argument("--transabyss-kmer", type=int, default=32, help="Trans-ABySS k-mer (default: 32)")
    p.add_argument(
        "--max-parallel", type=int, default=2,
        help="max concurrent jobs within each independent stage that benefits from "
             "it (orthofuser vs. merge/orthotransrate; transrate vs. strandeval), "
             "splitting --cpu/--mem across however many run at once; the 4 "
             "assemblers instead run as two fixed resource lanes (see "
             "TRINITY_LANE_SHARE), unaffected by this flag; other CPU-bound stages "
             "that don't benefit (diamond, orp_diamond, salmon, busco) always run "
             "sequentially at full --cpu; 1 disables concurrency for the remaining "
             "stages entirely (default: 2)",
    )
    p.add_argument("--dir", default=None, help="working directory (default: current directory)")
    return p.parse_args()


def main():
    args = parse_args()
    pipeline = Pipeline(args)
    try:
        pipeline.main()
    except subprocess.CalledProcessError as e:
        sys.exit(f"\n*** step failed: {' '.join(str(c) for c in e.cmd)} (exit {e.returncode}) ***")
    except FileNotFoundError as e:
        sys.exit(f"\n*** required command not found: {e.filename} ***")


if __name__ == "__main__":
    main()
