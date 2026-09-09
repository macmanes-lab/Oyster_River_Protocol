#!/usr/bin/env python3
"""chowder.py - bring your own assemblies, and let the ORP merge them.

Usage:
    chowder.py --assemblies trinity.fasta spades.fasta [more.fasta ...] \\
               --read1 R1.fq.gz --read2 R2.fq.gz --cpu 24 --mem 110 \\
               --runout myrun

Everything oyster.py does after its four assemblers have finished, run over
assemblies you already have: OrthoFinder groups the contigs, pytransrate
scores them, the best-scoring member of each orthogroup is kept, the
contigs that no orthogroup covered are rescued by their diamond hits,
cd-hit-est collapses what is left, and the result is quantified, filtered,
BUSCO'd, scored and reported exactly as a full ORP run would be.

It is the same code, not a copy of it: chowder subclasses oyster.py's
Pipeline and reuses `merge_and_report` wholesale. The selection rule, the
group ordering that reaches cd-hit-est's tie-breaks, and the scoring are
therefore the same ones by construction, and cannot drift from ORP's.

Two things to know before running it:

  * Reads are not optional. The merge scores every contig against them
    (pytransrate), quantifies the survivors (salmon), and strand-checks the
    result -- so a merge needs the library the assemblies came from, and
    trims and error-corrects it the way ORP always does. Pass
    --reads-are-corrected if you are handing over reads that have already
    been through trimmomatic and rcorrector.

  * Contig names are prefixed with the label of the assembly they came
    from. Two assemblies of the same library routinely share contig names
    -- two Trinity runs both start at TRINITY_DN0_c0_g1_i1 -- and every
    stage downstream is keyed by contig name, so the prefix is what keeps
    them apart. It doubles as provenance: every contig in the final
    assembly says which input it survived from.

Assembly order is significant, and it is the order you list them in:

  * it sets contig order in the pooled fasta, which reaches cd-hit-est,
    where input order breaks length ties -- so it can decide which of two
    equally long, 98%-identical contigs survives; and
  * build_list5.py keeps the *first* diamond hit per gene in that order, so
    for contigs no orthogroup covered, earlier assemblies are preferred.

List the assembly you trust most first.
"""

import argparse
import copy
import gzip
import re
import subprocess
import sys
import time
from pathlib import Path

from oyster import ASSEMBLER_TOOLS, CHECK_TOOLS, RED, RESET, Assembly, Pipeline

HERE = Path(__file__).resolve().parent

# Labels that would collide with a file the pipeline writes itself under
# assemblies/<runout>.*.
RESERVED_LABELS = frozenset({"orthomerged", "orp", "orp.intermediate", "filter", "flagstat"})


def sanitise_label(raw):
    """Turn a filename into a label safe in both filenames and contig IDs."""
    name = raw
    for suffix in (".gz", ".bz2"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    name = re.sub(r"\.(fa|fasta|fna|fas|seq)$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
    return name or "assembly"


def derive_labels(paths, explicit=None):
    """One unique, safe label per input assembly.

    Derived from the filename unless --labels was given. Duplicates get a
    numeric suffix rather than being rejected: two files called
    trinity.fasta in different directories is an ordinary way to hold two
    assemblies of the same library.
    """
    if explicit:
        if len(explicit) != len(paths):
            sys.exit(f"--labels: got {len(explicit)} for {len(paths)} assemblies")
        labels = [sanitise_label(l) for l in explicit]
    else:
        labels = [sanitise_label(Path(p).name) for p in paths]

    seen, out = {}, []
    for label in labels:
        if label.lower() in RESERVED_LABELS:
            label = f"{label}_input"
        n = seen.get(label.lower(), 0) + 1
        seen[label.lower()] = n
        out.append(label if n == 1 else f"{label}_{n}")
    return out


def open_maybe_gzip(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


class Chowder(Pipeline):
    """oyster.py's Pipeline, fed assemblies instead of building them."""

    def __init__(self, args):
        self.sources = [Path(p).resolve() for p in args.assemblies]
        labels = derive_labels(self.sources, args.labels)
        # One label per assembly, used for its file under assemblies/, its
        # diamond output, its unique-gene count, its line in the quality
        # report and the prefix on its contig names -- so a contig can be
        # traced from the final assembly back to the file it came from by
        # reading its name.
        # On a copy: Pipeline reads the assembly set off args, and quietly
        # rewriting the caller's namespace from paths to Assembly records
        # would make constructing a second pipeline from it fail.
        args = copy.copy(args)
        args.assemblies = [
            Assembly(f"{label}.fasta", label, label, label.upper()) for label in labels
        ]
        self.reads_are_corrected = args.reads_are_corrected
        super().__init__(args)
        self.ingest_done = self.assemblies_dir / f"{self.runout}.ingest.done"

    # -- preflight ---------------------------------------------------------

    def required_tools(self):
        """Everything a full run checks, less the assemblers, plus bwa.

        The three assemblers are the one part of the ORP chowder never
        reaches, so demanding them would fail a preflight over software this
        run has no use for. bwa takes their place on the list because
        strandeval maps with it out of the orp_trinity env, and that env is
        no longer proven present by the Trinity check.
        """
        tools = [t for t in CHECK_TOOLS if t not in ASSEMBLER_TOOLS]
        return tools + [("orp_trinity", "bwa", "BWA")]

    def readcheck(self):
        """Existence only -- there is no assembler k-mer to be too long for."""
        for label, path in (("READ1", self.read1), ("READ2", self.read2)):
            if not path.exists():
                sys.exit(f"\n\n\n\n ERROR: YOUR {label} FILE DOES NOT EXIST AT THE "
                         "LOCATION YOU SPECIFIED\n\n\n\n ")

    def assemblycheck(self):
        for src in self.sources:
            if not src.exists():
                sys.exit(f"\n*** assembly not found: {src} ***")
            if src.stat().st_size == 0:
                sys.exit(f"\n*** assembly is empty: {src} ***")
        if len(self.sources) < 2:
            sys.exit("\n*** chowder merges assemblies: give it at least two ***")

    def run_inputs(self):
        return [self.read1, self.read2] + self.sources

    def welcome(self):
        print(RED)
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("                    _.-~~~~~-._")
        print("                 .-~   o   o   ~-.")
        print("                (   .-'~~~~~'-.   )        OYSTER RIVER CHOWDER")
        print(f"                 '-.___________.-'         version {self.version}")
        print("                     '-.___.-'")
        print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + RESET + "\n")
        print(f"    Merging {len(self.assemblies)} assemblies, in this order:\n")
        for a, src in zip(self.assemblies, self.sources):
            print(f"      {a.diamond_label:<24} {src}")
        print("\n    Order matters: it breaks cd-hit-est's length ties and ranks the")
        print("    diamond rescue. Contigs are renamed <label>_<original>.\n")

    # -- ingest ------------------------------------------------------------

    def ingest(self):
        """Copy each input under assemblies/, prefixing every contig name.

        Prefixing is not cosmetic. Contig names are the key every downstream
        stage joins on -- OrthoFinder's orthogroups, pytransrate's
        contigs.csv, filter.py's keep-lists -- and two assemblies of one
        library routinely share them. Unprefixed, two contigs called
        TRINITY_DN0_c0_g1_i1 would be silently treated as one.
        """
        for a, src in zip(self.assemblies, self.sources):
            dst = self.assembly_fasta(a)
            prefix = f"{a.diamond_label}_"
            seen, n = set(), 0
            with open_maybe_gzip(src) as inf, open(dst, "w") as out:
                for line in inf:
                    if not line.startswith(">"):
                        out.write(line)
                        continue
                    n += 1
                    name = line[1:].split(None, 1)[0] if line[1:].strip() else f"contig{n}"
                    if name in seen:
                        sys.exit(
                            f"\n*** {src} contains the contig name {name!r} more than "
                            "once. Contig names have to be unique within an assembly: "
                            "every stage from OrthoFinder on joins on them. ***"
                        )
                    seen.add(name)
                    rest = line[1:][len(name):].rstrip("\n")
                    out.write(f">{prefix}{name}{rest}\n")
            if n == 0:
                sys.exit(f"\n*** no sequences found in {src} -- is it FASTA? ***")
            print(f"[ingest] {src} -> {dst}  ({n} contigs, renamed {prefix}*)")
            # Same deal as oyster.py's assembler lanes: start the gzip as
            # soon as the file stops being written, so cleanup() at the end
            # has nothing to do but unlink.
            self.compress_async(dst)
        self.ingest_done.write_text(
            "\n".join(f"{a.diamond_label}\t{src}" for a, src in zip(self.assemblies, self.sources)) + "\n"
        )

    # -- orchestration -----------------------------------------------------

    def main(self):
        pipeline_start = time.time()
        self.setup()
        if self.already_complete():
            return
        self.timing_init()
        self.check()
        self.welcome()
        self.readcheck()
        self.assemblycheck()
        self.step(
            "ingest", self.assembly_fasta_paths() + [self.ingest_done],
            self.sources, self.ingest, timed=False,
        )
        self.prepare_reads()
        self.merge_and_report(pipeline_start)

    def prepare_reads(self):
        """Trim and correct, unless told the reads are already corrected."""
        if not self.reads_are_corrected:
            return super().prepare_reads()
        c1, c2 = self.cor1(), self.cor2()
        self.rcorr_dir.mkdir(parents=True, exist_ok=True)
        for src, dst in ((self.read1, c1), (self.read2, c2)):
            if dst.exists() and not self.needs_run([dst], [src]):
                continue
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src.resolve())
        print(f"[reads] --reads-are-corrected: using {self.read1} / {self.read2} as they are")


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge assemblies you already have, using the ORP.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Assemblies are merged in the order given; list the one you trust most first.",
    )
    version = (HERE / "version.txt").read_text().strip()
    p.add_argument("--version", action="version", version=f"Oyster River Protocol {version} (chowder)")
    p.add_argument("--assemblies", nargs="+", required=True, metavar="FASTA",
                   help="two or more assemblies to merge, best first (.gz is fine)")
    p.add_argument("--labels", nargs="+", default=None, metavar="LABEL",
                   help="one name per assembly, used for its contig-name prefix and "
                        "its line in the quality report (default: from the filenames)")
    p.add_argument("--read1", required=True, help="path to R1 fastq(.gz)")
    p.add_argument("--read2", required=True, help="path to R2 fastq(.gz)")
    p.add_argument("--reads-are-corrected", action="store_true",
                   help="the reads have already been through trimmomatic and "
                        "rcorrector; skip both (default: off)")
    p.add_argument("--mem", type=int, default=110, help="memory in GB (default: 110)")
    p.add_argument("--cpu", type=int, default=16, help="CPU threads (default: 16)")
    p.add_argument("--busco-threads", type=int, default=None, help="BUSCO threads (default: same as --cpu)")
    p.add_argument("--runout", default="USER_RUN", help="run name prefix (default: USER_RUN)")
    p.add_argument("--lineage", default="eukaryota_odb12.2", help="BUSCO lineage (default: eukaryota_odb12.2)")
    p.add_argument("--tpm-filt", type=float, default=0, help="TPM filter threshold (default: 0)")
    p.add_argument(
        "--max-parallel", type=int, default=2,
        help="max concurrent jobs within each independent stage that benefits "
             "from it (orthofuser vs. merge/orthotransrate; transrate vs. "
             "strandeval), splitting --cpu/--mem across however many run at "
             "once; 1 disables it (default: 2)",
    )
    p.add_argument(
        "--keep-intermediates", action="store_true",
        help="keep every file the run produces: skips the end-of-run cleanup and "
             "the reclaim of the trimmed reads. For debugging (default: off)",
    )
    p.add_argument("--dir", default=None, help="working directory (default: current directory)")
    return p.parse_args()


def main():
    args = parse_args()
    pipeline = Chowder(args)
    try:
        pipeline.main()
    except subprocess.CalledProcessError as e:
        sys.exit(f"\n*** step failed: {' '.join(str(c) for c in e.cmd)} (exit {e.returncode}) ***")
    except FileNotFoundError as e:
        sys.exit(f"\n*** required command not found: {e.filename} ***")
    finally:
        pipeline.finish_compression()


if __name__ == "__main__":
    main()
