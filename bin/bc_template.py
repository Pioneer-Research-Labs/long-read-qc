#!/usr/bin/env python

from Bio import SeqIO
import sys
import re
from snapgene_reader import snapgene_file_to_seqrecord

def extract_flanks(path, out_type, trim_insert_flank=False,  construct=None):
    """
    Extracts flanking sequences from a genbank file and returns them in the format required by cutadapt
    :param path: Path to the genbank file containing the flanking sequences
    :param out_type: string representing the annotation type to be extracted. Either 'cutadapt_barcode' or 'cutadapt_insert'
        :param trim_insert_flank: boolean indicating whether any trimming should occur for flanking insert sequences
    :param construct: string representing the construct name, used for bartender output
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
        seq_up = str(ordered_seqs[0].seq)
        seq_dn = str(ordered_seqs[1].seq)
        if trim_insert_flank:

            seq_up = trim_string(seq_up, 0, 35)
            seq_dn = trim_string(seq_dn, 35, 0)
        bc = f'{seq_up}...{seq_dn}'
    elif out_type == "cutadapt_genome_tag":
        ordered_seqs = order_based_on_position(filter_seqs_on_name(seqs, r'Barcode_(R|F)|INSERT(UP|DN)|Secondary_Barcode_for_Donor_gDNA'))
        if len(ordered_seqs) != 5:
            raise ValueError("Not enough sequences found for cutadapt_genome_tag. Expected at least 5 sequences, found {}".format(len(ordered_seqs)))
        # We're assuming when ordered, the 2nd and 3rd are the INSERTDN and Secondary Barcode for Donor gDNA sequences

        # ---<INSERTUP>---<INSERTDN---<Secondary Barcode for Donor gDNA>-->---<Barcode F>---<Barcode R>---

        # We care about the span between INSERTDN and Secondary Barcode for Donor gDNA and the span between
        # the end of the Secondary Barcode for Donor gDNA and the start of Barcode F
        #                 <--<INSERTDN>---<Secondary Barcode for Donor gDNA>-->---<Barcode F>

        # The first adapter is formed by truncating the insert dn sequence to the start of the gDNA tag sequence
        # The second adapter is formed by truncating the end of gDNA tag to the start of the Barcode F sequence

        insert_dn = ordered_seqs[1]
        g_tag = ordered_seqs[2]
        barcode_f = ordered_seqs[3]
        insert_up_start = insert_dn.features[0].location.start
        g_tag_start = g_tag.features[0].location.start
        g_tag_end = g_tag.features[0].location.end
        barcode_f_start = barcode_f.features[0].location.start
        adapter_1 = str(insert_dn.seq[:g_tag_start - insert_up_start])
        intervening_seq = get_intervening_sequence(construct, g_tag_end, barcode_f_start)

        adapter_2 = str(intervening_seq)

        bc = f'{adapter_1}...{adapter_2}'

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

def trim_seq(five_prime_trim_length, three_prime_trim_length):
    """

    :param five_prime_trim_length: int representing number of bases to trim from 5' end
    :param three_prime_trim_length: int representing number of bases to trim from 3' end
    :return: trimmed sequence as a string
    """

def trim_string(input_string, start_trim, end_trim):
    """
    Trims characters from the start and end of a string.

    :param input_string: The string to be trimmed
    :param start_trim: Number of characters to trim from the start
    :param end_trim: Number of characters to trim from the end
    :return: The trimmed string
    """
    return input_string[start_trim:len(input_string) - end_trim]


def get_intervening_sequence(construct, start, end):
    """
    Extracts the intervening sequence between two positions in a sequence.

    :param construct: The path to the SnapGene file
    :param start: Start position (0-based)
    :param end: End position (0-based)
    :return: The intervening sequence as a string
    """
    seq = snapgene_file_to_seqrecord(construct)
    return seq[start:end].seq

if __name__ == "__main__":

    if len(sys.argv) == 3:
        extract_flanks(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        extract_flanks(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 5:
        extract_flanks(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])