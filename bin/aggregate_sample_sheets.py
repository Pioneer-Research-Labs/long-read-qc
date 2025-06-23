#!/usr/bin/env python
import sys
import pandas as pd

def aggregate(sample_sheet_map):
    """
    Concatenate files listed in a sample sheet map into a single file.
    """
    # Read file paths from the sample sheet map
    with open(sample_sheet_map) as f:
        # Create a dictionary of sample names and file paths
        sample_file_dict = {line.split('\t')[0]: line.split('\t')[1].rstrip('\n') for line in f}

    df_to_concatenate = []

    for key,val in sample_file_dict.items():
        df = pd.read_csv(val, engine='c', header=None)
        df_to_concatenate.append(df)

    concatenated_df = pd.concat(df_to_concatenate, ignore_index=True)

    concatenated_df.columns =['id', 'genome','construct','file']

    # Write to output
    concatenated_df.to_csv('aggregated_sample_sheet.csv', index=False, header=True )

if __name__ == "__main__":
    aggregate(sys.argv[1])