#!/usr/bin/env python

from snapgene_reader import snapgene_file_to_dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import re


def extract_flanks(path):
    """
    Parse the SnapGene file and extract the features matching the given patterns.
    :param path: Path to the SnapGene file
    :return: Written Genbank file with the flanking sequences
    """
    con = snapgene_file_to_dict(path)
    part_labels = [r'INSERT(UP|DN)', r'BARCODE[0-9]{0,2}(UP|DN)',  r'Barcode (R|F)', 'Secondary Barcode for Donor gDNA']
    part_list = {}
    for x in con['features']:
        if any(re.match(reg, x['name']) for reg in part_labels):
            part_list[x['name']] = (x['start'], x['end'])

    flanking_seqs = []
    for x,y in part_list.items():
        seq = con['seq'][y[0]:y[1]]
        # We need location information later on to correctly orient adapters in cutadapt
        f = FeatureLocation(y[0], y[1])
        feature = SeqFeature(f, type = "misc_feature")
        # Locus names are not allowed to have spaces, so we replace them with underscores
        record = SeqRecord(Seq(seq), id = x.replace(" ","_"), description = "", features=[feature], annotations={"molecule_type": "DNA"})
        flanking_seqs.append(record)

    return flanking_seqs


if __name__ == "__main__":
    construct_file = sys.argv[1]
    flanks = extract_flanks(construct_file)
    # Write Genbank format to preserve location information
    SeqIO.write(flanks, "flanking.gb", "genbank")