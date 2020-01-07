#!/usr/bin/python

import sys
from Bio import SeqIO


#usage: python long.seq.py in.fasta out.fasta 200

INPUT_SEQ_ITERATOR = SeqIO.parse(sys.argv[1], "fasta")
SHORT_SEQ_ITERATOR = (record for record in INPUT_SEQ_ITERATOR \
                      if len(record.seq) > int(sys.argv[3]))

OUTPUT_HANDLE = open(sys.argv[2], "w")
SeqIO.write(SHORT_SEQ_ITERATOR, OUTPUT_HANDLE, "fasta")
OUTPUT_HANDLE.close()
