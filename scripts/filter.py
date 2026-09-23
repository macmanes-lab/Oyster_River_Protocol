"""
%prog some.fasta wanted-list.txt

Write the records of some.fasta whose ID is listed in wanted-list.txt to
stdout, in some.fasta's order.

This used to be Biopython's SeqIO with `seq.id in wanted` over a *list*, a
linear scan per record: millions of records in merged.fasta against every
line of good.<run>.list made it quadratic, and on big assemblies it was the
slow step. The IDs are a set now, and the parse is plain bytes so that the
records not wanted -- nearly all of them -- are skipped without being built.

The output is byte-identical to what SeqIO.parse/SeqIO.write produced, and
that is deliberate: this is orthomerged.fasta, whose record order and text
reach cd-hit-est. So, as Biopython does: the ID is the first word of the
header; the header is written back with trailing whitespace stripped; the
sequence has spaces, tabs and carriage returns removed and is re-wrapped at
60 columns; a record with no sequence is written as its header alone; and a
record listed once but present twice is written twice.
"""
import sys

WRAP = 60


def filter_fasta(fasta, wanted, out):
    keep = False
    seq = []

    def flush():
        data = b"".join(seq).translate(None, b" \t\r\n")
        for i in range(0, len(data), WRAP):
            out.write(data[i:i + WRAP])
            out.write(b"\n")

    for line in fasta:
        if line[:1] == b">":
            if keep:
                flush()
            title = line[1:].rstrip()
            words = title.split(None, 1)
            keep = (words[0] if words else b"") in wanted
            if keep:
                out.write(b">" + title + b"\n")
                seq = []
        elif keep:
            seq.append(line)
    if keep:
        flush()


def main():
    with open(sys.argv[2], "rb") as fh:
        wanted = {line.strip() for line in fh}
    with open(sys.argv[1], "rb") as fasta:
        first = fasta.readline()
        if first and first[:1] != b">":
            sys.exit(f"{sys.argv[1]}: does not start with a '>' header line")
        fasta.seek(0)
        filter_fasta(fasta, wanted, sys.stdout.buffer)


if __name__ == "__main__":
    main()
