#!/usr/bin/env python3

#usage: python strandeval.py --assembly test.fasta --read1 R1.fq --read2 R2.fq
#           [--runout test] [--cpu 16] [--pairs 400000] [--dir .] [--no-cleanup]
#
#The Oyster River Strand Exam Tool on its own, for any assembly and read pair:
#the strandeval step oyster.py and chowder.py run on <run>.ORP.fasta, without
#the run around it. It replaces strandeval.mk.
#
#Samples --pairs read pairs (seqtk, fixed seed), maps them to the assembly
#(bwa mem), counts the strand each transcript's first reads land on
#(examine_strand.pl), and prints the histogram of (plus - minus) / total
#across transcripts. docs/strandexamine.md says how to read it.
#
#The work is oyster.py's strand_exam(), not a copy of it, so this gives the
#same answer the pipeline would. Like the pipeline, it needs bwa, seqtk, hist
#and Trinity in the orp_trinity conda env and samtools in orp.
#
#Writes, under --dir:
#  reports/<runout>.strandeval_summary.txt   the histogram, as printed
#  reports/<runout>.flagstat                 samtools flagstat of the sample
#The BAM, bwa index and per-transcript table (<runout>.dat) are deleted when
#it finishes, unless --no-cleanup.

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from oyster import STRAND_SAMPLE_PAIRS, line_buffer_stdio, strand_exam  # noqa: E402

TOOLS = (("orp_trinity", ("bwa", "seqtk", "hist", "Trinity")),
         ("orp", ("samtools",)))


def run(cmd, cwd=None):
    print("+", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def missing_tools():
    """'env/tool' for each tool strand_exam needs that its env lacks."""
    missing = []
    for env, binaries in TOOLS:
        try:
            result = subprocess.run(
                ["conda", "run", "-n", env, "which", *binaries],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
        except FileNotFoundError:
            return ["conda"]
        found = {Path(line.strip()).name for line in result.stdout.splitlines()}
        missing += [f"{env}/{b}" for b in binaries if b not in found]
    return missing


def main():
    p = argparse.ArgumentParser(description="Oyster River Strand Exam Tool: how stranded "
                                            "are the reads, as mapped to this assembly?")
    p.add_argument("--assembly", required=True, help="assembly fasta")
    p.add_argument("--read1", required=True, help="R1 fastq(.gz)")
    p.add_argument("--read2", required=True, help="R2 fastq(.gz)")
    p.add_argument("--runout", default=None,
                   help="prefix for the files written (default: the assembly's name)")
    p.add_argument("--cpu", type=int, default=16, help="threads (default: 16)")
    p.add_argument("--pairs", type=int, default=STRAND_SAMPLE_PAIRS,
                   help=f"read pairs to sample (default: {STRAND_SAMPLE_PAIRS}, as the pipeline)")
    p.add_argument("--dir", default=None, help="working directory (default: current directory)")
    p.add_argument("--no-cleanup", action="store_true",
                   help="keep the BAM, bwa index and per-transcript table (<runout>.dat)")
    args = p.parse_args()
    line_buffer_stdio()

    assembly, r1, r2 = (Path(f).resolve() for f in (args.assembly, args.read1, args.read2))
    for f in (assembly, r1, r2):
        if not f.is_file():
            sys.exit(f"*** no such file: {f} ***")
    runout = args.runout or assembly.name.split(".")[0]
    workdir = Path(args.dir).resolve() if args.dir else Path.cwd()
    reports = workdir / "reports"
    reports.mkdir(parents=True, exist_ok=True)

    missing = missing_tools()
    if missing:
        sys.exit(f"*** not installed: {', '.join(missing)} ***")

    try:
        summary = strand_exam(run, assembly, r1, r2, runout, workdir, args.cpu,
                              flagstat_out=reports / f"{runout}.flagstat",
                              pairs=args.pairs, keep_scratch=args.no_cleanup)
    except subprocess.CalledProcessError as e:
        sys.exit(f"\n*** step failed: {' '.join(str(c) for c in e.cmd)} (exit {e.returncode}) ***")
    print(summary)
    (reports / f"{runout}.strandeval_summary.txt").write_text(summary)


if __name__ == "__main__":
    main()
