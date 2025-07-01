#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pathlib import Path

def plot(inserts_full, insert_truncated, all_seq_lengths, length_untrimmed, sample):
    """
    Plot comparison of full and truncated inserts, along with read lengths.
    Args: inserts_full (str): Path to the file with full insert and lengths.
          insert_truncated (str): Path to the file with truncated inserts and lengths.
          all_seq_lengths (str): Path to the file with all sequence lengths for sample.
          length_untrimmed (str): Path to the file with lengths of untrimmed reads.
          sample (str): Sample name for labeling the plots.
    Returns: None
    """
    # Read the full and truncated inserts
    full_df = pd.read_csv(inserts_full, sep='\t', names=['id', 'length'])
    full_df['source'] = 'intact_flanks'
    truncated_df = pd.read_csv(insert_truncated, sep='\t', names=['id', 'length'])
    truncated_df['source'] = 'truncated_flanks'

    # If empty dataframes are returned, return early
    if full_df.empty or truncated_df.empty:
        print("No data to plot. Exiting.")
        plt.savefig('truncated_vs_intact_flanks_comparison_plot.png')
        return
    all_inserts_df = pd.concat([truncated_df, full_df])

    # Read lengths for intact and truncated inserts
    all_seq_lengths_df = pd.read_csv(all_seq_lengths, sep='\t', names=['id', 'read_length'])
    all_inserts_with_lengths = pd.merge(all_inserts_df, all_seq_lengths_df, on='id', how='left')
    combined_with_read_length = all_inserts_with_lengths[all_inserts_with_lengths['read_length'].notna()]
    no_inserts = pd.read_csv(length_untrimmed, sep='\t', names=['id', 'read_length'])


    fig, ax = plt.subplots(ncols=4, figsize=(25, 6))
    sns.histplot(all_inserts_df,
                 x='length',
                 hue='source',
                 ax=ax[0])
    ax[0].set_xlabel('Insert length')
    ax[0].set_title(f'Insert lengths for {sample}')
    sns.countplot(x="source", data=all_inserts_df, dodge=False, ax=ax[1])
    ax[1].set_title(f'Raw counts of inserts for {sample}')

    sns.histplot(combined_with_read_length,
                 x='read_length',
                 hue='source',
                 ax=ax[2])
    ax[2].set_xlabel('Read length')
    ax[2].set_title(f'Read lengths for trimmed seqs in {sample}')
    # Plot read lengths for reads with no inserts

    sns.histplot(no_inserts,
                 x='read_length',
                 ax=ax[3])
    ax[3].set_xlabel('Read length')
    ax[3].set_title(f'Read lengths for untrimmed reads from {sample}')
    plt.tight_layout()

    output_file = Path('truncated_vs_intact_flanks_comparison_plot.png')
    plt.savefig(output_file)
    plt.close()

if __name__ == '__main__':
    full_insert_lengths = sys.argv[1]
    truncated_insert_lengths = sys.argv[2]
    all_sequence_lengths = sys.argv[3]
    untrimmed_insert_lengths = sys.argv[4]
    sample_name = sys.argv[5]
    plot( full_insert_lengths, truncated_insert_lengths, all_sequence_lengths, untrimmed_insert_lengths, sample_name)
