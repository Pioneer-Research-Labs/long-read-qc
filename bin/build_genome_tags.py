#!/usr/bin/env python
import sys
import pandas as pd
from snapgene_reader import snapgene_file_to_dict
import re
from Bio.Seq import Seq

def extract_upstream_bases(snapgene_file, bp_prefix_length=5):
    """
    Extracts the sequence upstream of the 'Secondary Barcode for Donor gDNA' feature
    in a SnapGene file.
    :param snapgene_file: Path to the SnapGene file
    :param bp_prefix_length: Number of base pairs to include upstream of the feature
    :return: The upstream sequence as a string
    """
    construct = snapgene_file_to_dict(snapgene_file)
    tag = ['Secondary Barcode for Donor gDNA']
    upstream_sequence = None
    for feature in construct['features']:
        if any(re.match(reg, feature['name']) for reg in tag):
            # Extract the sequence upstream of the feature
            start = feature['start'] - bp_prefix_length
            upstream_sequence = construct['seq'][start:feature['start']]
    return upstream_sequence


def build_genome_tags(genome_tag_file, upstream_bases):
    """

    """
    tesseract_df = pd.read_csv(genome_tag_file)
    tesseract_df['8bp Barcode'] = tesseract_df['8bp Barcode'].str.upper()
    # Prefix the 8bp barcode with the upstream bases
    tesseract_df['Base Strain'] = tesseract_df['Base Strain'].str.replace(" ", "_")
    # Create the reverse complement of the 8bp barcode
    tesseract_df['8bp Barcode RC'] = tesseract_df['8bp Barcode'].apply(lambda x: str(Seq(x).reverse_complement()))
    # Write the 8bp barcode, its reverse complement, and the Base Strain to a new file
    columns_to_keep = ['8bp Barcode', '8bp Barcode RC', 'Base Strain']
    tesseract_df.to_csv("prefixed_tags.csv", columns=columns_to_keep)
    # Combine barcode and it's reverse complement into one column
    combined_df = pd.concat([tesseract_df['8bp Barcode'],tesseract_df['8bp Barcode RC']], axis=0).reset_index(drop=True)
    combined_df.to_csv("combined_tags.csv", header=False, index=False)

if __name__ == "__main__":
    sequence = extract_upstream_bases(sys.argv[1])
    build_genome_tags(sys.argv[2], sequence)