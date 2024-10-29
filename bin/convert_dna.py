#!/usr/bin/env python
from snapgene_reader import snapgene_file_to_seqrecord
from Bio import SeqIO
import sys
import os


def convert_dna(path):
    seq = snapgene_file_to_seqrecord(path)
    base = os.path.splitext(os.path.basename(path))[0]
    seq.id = base
    seq.description = ""

    return seq


seq = convert_dna(sys.argv[1])
SeqIO.write(seq, sys.stdout, 'fasta-2line')