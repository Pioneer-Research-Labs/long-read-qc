#!/usr/bin/env python

from Bio import SeqIO
import sys
import re

def extract_flanks(path, out_type):
    """
    Extracts flanking sequences from a genbank file and returns them in the format required by cutadapt
    :param path: Path to the genbank file containing the flanking sequences
    :param out_type: string representing the annotation type to be extracted. Either 'cutadapt_barcode' or 'cutadapt_insert'
    :return: string representing the flanking sequences in the format required by cutadapt
    """

    if out_type == "bartender":
        mid = sys.argv[3]

    seqs = [x for x in SeqIO.parse(path, 'genbank')]

    if out_type == "cutadapt_barcode":
        ordered_seqs =  order_based_on_position(filter_seqs_on_name(seqs, r'BARCODE[0-9]{0,2}(UP|DN)'))
        bc = f'{str(ordered_seqs[0].seq)}...{str(ordered_seqs[1].seq)}'
    elif out_type == "cutadapt_insert":
        ordered_seqs =  order_based_on_position(filter_seqs_on_name(seqs, r'INSERT(UP|DN)'))
        bc = f'{str(ordered_seqs[0].seq)}...{str(ordered_seqs[1].seq)}'
    elif out_type == "bartender":
        ordered_seqs = order_based_on_position(filter_seqs_on_name(seqs, r'BARCODE[0-9]{0,2}(UP|DN)'))
        bc = f'{str(ordered_seqs[0].seq[-5:])}[20]{mid}[20]{str(ordered_seqs[1].seq[:5])}'
    else:
        raise ValueError("Please provide a type of either 'cutadapt' or 'bartender'")

    sys.stdout.write(bc)

def filter_seqs_on_name(seqs, pattern):
    """
    Filters the sequences based on the name
    :param seqs: List of SeqRecords
    :param pattern: Regex name of the sequence
    :return: List of SeqRecords with the given name
    """
    return [x for x in seqs if  re.match(pattern, x.name)]

def order_based_on_position(seqs):
    """
    Orders the sequences based on their position
    :param seqs: List of SeqRecords
    :return: List of SeqRecords ordered by position
    """
    return sorted(seqs, key=lambda x: x.features[0].location.start)

if __name__ == "__main__":
    extract_flanks(sys.argv[1], sys.argv[2])