#!/usr/bin/env python

import pandas as pd
import sys

from bin.combine_insert_coverage import write_empty_file


def concatenate_barcode_count_files(file_map):
    with open(file_map, 'r') as f:
        sample_file_dict = {line.split()[0]: line.split()[1] for line in f}
    df_to_concatenate = []
    for key,val in sample_file_dict.items():
        try:
            counts = pd.read_table(val, names=['barcode_seq', 'barcode_count']).assign(sample=key)
            df_to_concatenate.append(counts)
        except pd.errors.EmptyDataError:
            print(f'No barcode data for {key} found in {val}')
            continue

    if not df_to_concatenate:
        write_empty_file('concatenated_barcode_counts.csv')
        return None
    sample_df = pd.concat(df_to_concatenate)
    sample_df.to_csv('concatenated_barcode_counts.csv', index=False)
    return sample_df


if __name__ == '__main__':
    barcode_file_map = sys.argv[1]
    barcode_counts = concatenate_barcode_count_files(barcode_file_map)



