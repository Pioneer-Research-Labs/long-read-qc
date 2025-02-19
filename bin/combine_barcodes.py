#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys


def write_empty_file(filename):
    with open(filename, "w") as my_empty_csv:
        pass


def plot_barcode_length_histogram(barcode_df):
    df_plot = barcode_df
    fig, ax = plt.subplots(figsize=(10, 6))
    if df_plot is None:
        plt.savefig('barcode_length_distribution.png')  # empty plot
        write_empty_file('barcode_lengths.csv')
        return
    df_plot.to_csv('barcode_lengths.csv', index=False)

    sns.boxplot(data=df_plot, x="sample", y="barcode_len", ax=ax)
    ax.set_title('Distribution of barcode lengths')
    plt.savefig('barcode_length_distribution.png', bbox_inches='tight')


def plot_proportions_of_barcodes(barcode_df):
    if barcode_df is None:
        plt.savefig('barcode_proportions.png')
        write_empty_file('barcode_proportions.csv')
        return
    # Categorize barcodes by length and plot fraction composition for each library
    # Use 1bp tolerance on either end
    barcode_length = 44  # expected bp length of barcodes
    empty_cutoff = 5  # cutoff for labeling a vector as empty
    barcode_df['barcode_exists'] = 'other'
    barcode_df.loc[
        barcode_df.barcode_len.between(barcode_length - 1, barcode_length + 1), 'barcode_exists'] = 'true'
    barcode_df.loc[barcode_df.barcode_len <= empty_cutoff, 'barcode_exists'] = 'false'

    df_plot = barcode_df.groupby('sample')['barcode_exists'].value_counts(normalize=True).to_frame()
    barcode_df.to_csv('barcode_proportions.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.barplot(df_plot,
                y='sample',
                x='proportion',
                hue='barcode_exists',
                ax=ax)
    ax.set_title('Proportion of true/empty/other detected barcodes')
    plt.savefig('barcode_proportions.png', bbox_inches='tight')


def plot_copy_number(barcode_df):
    if barcode_df is None:
        plt.savefig('barcode_copy_number.png')  # empty plot
        write_empty_file('barcode_copy_number.csv')
        return
    df_plot = barcode_df.groupby('sample')['barcode_seq'].value_counts().to_frame().rename(
        columns={'count': 'counts_per_barcode'})

    fig, ax = plt.subplots(figsize=(10, 6))

    df_plot.to_csv('barcode_copy_number.csv', index=False)
    sns.histplot(df_plot,
                 x='counts_per_barcode',
                 hue='sample',
                 ax=ax)
    ax.set_title('Copy number of each unique barcode')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig('barcode_copy_number.png', bbox_inches='tight')


def concatenate_barcode_files(file_map):
    with open(file_map, 'r') as f:
        sample_file_dict = {line.split()[0]: line.split()[1] for line in f}
    df_to_concatenate = []
    for key,val in sample_file_dict.items():
        try:
            barcode_df = pd.read_table(val, names=['read', 'barcode_seq', 'barcode_len'], usecols=[0, 1, 3])
            barcode_df['sample'] = key
            df_to_concatenate.append(barcode_df)
        except pd.errors.EmptyDataError:
            print(f'No barcode data for {key} found in {val}')
            continue

    if not df_to_concatenate:
        write_empty_file('concatenated_barcodes.csv')
        return None
    barcodes_df = pd.concat(df_to_concatenate)
    barcodes_df.to_csv('concatenated_barcodes.csv', index=False)
    return barcodes_df


if __name__ == '__main__':
    # Load barcodes from all samples
    barcodes = concatenate_barcode_files(sys.argv[1])
    # Plot barcode length histogram
    plot_barcode_length_histogram(barcodes)
    # Plot proportions of barcodes
    plot_proportions_of_barcodes(barcodes)
    # Plot copy number
    plot_copy_number(barcodes)