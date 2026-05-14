#!/usr/bin/env python

import pandas as pd

def empty_insert_plots(barcode_file, insert_file, title):
    barcodes = pd.read_csv(barcode_file, sep='\t', names=['id', 'seq', 'quality', 'length'])
    inserts = pd.read_csv(insert_file, sep='\t', names=['id', 'seq', 'quality', 'length'])

    with open("counts_of_inserts_barcodes_flanking_site_sequences.csv", "w") as summary_file:
        summary_file.write(
            f"Sample name,"
            f"Count of sequences with inserts,"
            f"Count of sequences with barcodes\n"
        )
        summary_file.write(
            f"{title},"
            f"{len(inserts)},"
            f"{len(barcodes)}\n"
        )

if __name__ == '__main__':
    import sys
    empty_insert_plots(sys.argv[1], sys.argv[2], sys.argv[3])
