#!/usr/bin/env python

from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import re

def extract_flanks(path):
    con = snapgene_file_to_dict(path)
    part_labels = [r'INSERT(UP|DN)', r'BARCODE[0-9]{0,2}(UP|DN)']
    part_list = {}
    for x in con['features']:
        if any(re.match(reg, x['name']) for reg in part_labels):
            part_list[x['name']] = (x['start'], x['end'])

    flanking_seqs = []
    for x,y in part_list.items():
        seq = con['seq'][y[0]:y[1]]
        flanking_seqs.append(SeqRecord(Seq(seq), id = x, description = ""))

    return flanking_seqs



construct_file = sys.argv[1]
flanks = extract_flanks(construct_file)
SeqIO.write(flanks, "flanking.fasta", "fasta-2line")