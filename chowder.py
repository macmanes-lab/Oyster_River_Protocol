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
    --corrected-reads if you are handing over reads that have already
    been through trimmomatic and rcorrector.

  * Contig names are prefixed with the label of the assembly they came
    from. Two assemblies of the same library routinely share contig names
    -- two Trinity runs both start at TRINITY_DN0_c0_g1_i1 -- and every
    stage downstream is keyed by contig name, so the prefix is what keeps
    them apart. It doubles as provenance: every contig in the final
    assembly says which input it survived from.

Assembly order is significant, and by default it is not yours to get wrong.
It matters twice:

  * it sets contig order in the pooled fasta, which reaches cd-hit-est,
    where input order breaks length ties -- so it can decide which of two
    equally long, 98%-identical contigs survives; and
  * build_list5.py keeps the *first* diamond hit per gene in that order, so
    for contigs no orthogroup covered, earlier assemblies are preferred.

Rather than let the sequence you happened to type decide either of those,
the assemblies are sorted by label and then permuted with a fixed seed, so
the order depends on the *set* of assemblies and not on how they were
listed. It is a seeded shuffle and not a random one on purpose: randomising
outright would mean the same command gave a different assembly on a
different day, and a resumed run disagreeing with the run it resumed. The
order used, and the seed that produced it, are printed at the start of the
run and written to assemblies/<run>.ingest.done.

This makes the order arbitrary and reproducible. It does not make the
pipeline order-independent -- the picks still depend on the order, and
being genuinely independent of it would mean breaking cd-hit-est's ties and
the rescue ranking on merit rather than on position. To measure what the
order is worth on your own data, run it twice with different --seed values
and diff the assemblies. To rank the assemblies yourself -- best first --
pass --assembly-order given.
"""

import argparse
import copy
import gzip
import random
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

# The seed strandeval already samples reads with, reused rather than adding a
# second arbitrary constant to the repo.
DEFAULT_SEED = 23894


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

    labels = [f"{l}_input" if l.lower() in RESERVED_LABELS else l for l in labels]

    # Two inputs can land on the same label -- two files called
    # trinity.fasta in different directories is an ordinary way to keep two
    # assemblies of one library. Numbering them in the order they were typed
    # would put the label assignment back under the user's control, and with
    # it the order (which is sorted by label), so the group is numbered by
    # source path instead: same set of files, same labels, however they were
    # listed. Every member of a colliding group is suffixed, including the
    # first -- a bare `trinity` beside a `trinity_2` reads as though the
    # bare one were somehow the real one.
    counts = {}
    for label in labels:
        counts[label.lower()] = counts.get(label.lower(), 0) + 1

    groups = {}
    for i, label in enumerate(labels):
        groups.setdefault(label.lower(), []).append(i)

    out = list(labels)
    for key, members in groups.items():
        if counts[key] == 1:
            continue
        for n, i in enumerate(sorted(members, key=lambda i: str(paths[i])), start=1):
            out[i] = f"{labels[i]}_{n}"
    return out


def shuffled_order(pairs, seed):
    """Put the assemblies in an order that doesn't depend on how they were typed.

    Order matters twice over downstream -- it breaks cd-hit-est's length
    ties and it ranks the diamond rescue -- and for assemblies nobody has
    ranked against each other, having the command line's order silently
    decide that is worse than having nothing decide it. So the order is
    taken out of the user's hands: sort by label first, so the result
    depends on the *set* of assemblies and not on the sequence they were
    listed in, then permute with a fixed seed.

    Deliberately a seeded shuffle and not a random one. Randomising outright
    would mean the same command produced a different assembly on a different
    day, and a resumed run disagreeing with the run it resumed -- trading a
    decision nobody made for one nobody can reproduce. With the seed fixed
    and recorded, the order is arbitrary but stable, and `--seed` makes the
    dependence measurable: run it twice with different seeds and the
    difference between the two assemblies is the size of the effect.

    Note what this does *not* do: the picks still depend on the order. It
    stops that order being an accident of typing; it does not make the
    pipeline order-independent, which would mean breaking cd-hit-est's ties
    and the rescue ranking on merit rather than on position.
    """
    ordered = sorted(pairs, key=lambda pair: pair[0])
    random.Random(seed).shuffle(ordered)
    return ordered


def open_maybe_gzip(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


class Chowder(Pipeline):
    """oyster.py's Pipeline, fed assemblies instead of building them."""

    def __init__(self, args):
        sources = [Path(p).resolve() for p in args.assemblies]
        labels = derive_labels(sources, args.labels)
        # One label per assembly, used for its file under assemblies/, its
        # diamond output, its unique-gene count, its line in the quality
        # report and the prefix on its contig names -- so a contig can be
        # traced from the final assembly back to the file it came from by
        # reading its name.
        self.seed = args.seed
        self.order_mode = args.assembly_order
        pairs = shuffled_order(list(zip(labels, sources)), self.seed) \
            if self.order_mode == "shuffled" else list(zip(labels, sources))
        self.sources = [src for _, src in pairs]
        # On a copy: Pipeline reads the assembly set off args, and quietly
        # rewriting the caller's namespace from paths to Assembly records
        # would make constructing a second pipeline from it fail.
        args = copy.copy(args)
        args.assemblies = [
            Assembly(f"{label}.fasta", label, label, label.upper()) for label, _ in pairs
        ]
        self.corrected_reads = args.corrected_reads
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
        if self.order_mode == "shuffled":
            origin = f"shuffled, seed {self.seed} -- not the order you listed them in"
        else:
            origin = "--assembly-order given: your order, so it ranks them"
        print(f"    Merging {len(self.assemblies)} assemblies ({origin}):\n")
        for a, src in zip(self.assemblies, self.sources):
            print(f"      {a.diamond_label:<24} {src}")
        print("\n    Order matters: it breaks cd-hit-est's length ties and ranks the")
        print("    diamond rescue, so it is recorded above and in reports/. Contigs")
        print("    are renamed <label>_<original>.\n")

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
        # The order is part of the result, not a detail of the invocation:
        # it decided cd-hit-est's ties and the rescue ranking, so a run has
        # to say which order it used and what produced it.
        header = (f"# assembly order: shuffled, seed {self.seed}"
                  if self.order_mode == "shuffled" else
                  "# assembly order: as given on the command line")
        self.ingest_done.write_text(
            header + "\n"
            + "\n".join(f"{i}\t{a.diamond_label}\t{src}"
                        for i, (a, src) in enumerate(zip(self.assemblies, self.sources), start=1))
            + "\n"
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
        if not self.corrected_reads:
            return super().prepare_reads()
        c1, c2 = self.cor1(), self.cor2()
        self.rcorr_dir.mkdir(parents=True, exist_ok=True)
        for src, dst in ((self.read1, c1), (self.read2, c2)):
            if dst.exists() and not self.needs_run([dst], [src]):
                continue
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src.resolve())
        print(f"[reads] --corrected-reads: using {self.read1} / {self.read2} as they are")


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
    p.add_argument("--assembly-order", choices=["shuffled", "given"], default="shuffled",
                   help="'shuffled' (default) puts the assemblies in a seeded, "
                        "reproducible order that does not depend on the order you "
                        "list them in; 'given' uses your order, which then ranks "
                        "them -- earlier assemblies win cd-hit-est's length ties "
                        "and are preferred by the diamond rescue")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                   help=f"seed for --assembly-order shuffled (default: {DEFAULT_SEED}). "
                        "Two runs differing only in this seed differ by exactly the "
                        "amount assembly order is worth")
    p.add_argument("--corrected-reads", action="store_true",
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
