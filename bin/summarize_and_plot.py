#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import csv
import plotly.io as pio
from pathlib import Path

# Color settings for seaborn and plotly
sns.set_theme(font_scale=.8)
sns.set_style('darkgrid')
pioneer_colors = ['#FF8633', '#423759', '#314942', '#FFA632', '#F7F3ED']
sns.set_palette(sns.color_palette(pioneer_colors))

pio.templates['pioneer'] = pio.templates["seaborn"]
pio.templates['pioneer'].layout.colorway = pioneer_colors
pio.templates.default = 'pioneer'


def write_empty_file(filename):
    """
    Write an empty file to the specified filename
    :param filename: str representing the filename
    """
    with open(filename, "w") as _:
        pass


def plot_full_genes_per_fragment(insert_cov_full):
    """
    Plot the number of full genes per insert
    :param insert_cov_full: DataFrame containing the number of full genes per insert

    """
    if insert_cov_full.empty:
        plt.savefig('full_genes_per_fragment.png') # empty plot
        write_empty_file('full_genes_per_fragment.csv')
        return
    insert_cov_full.to_csv('full_genes_per_fragment.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(insert_cov_full,
                 x='count',
                 hue='sample',
                 stat='frequency',
                 common_norm=True,
                 ax=ax)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_title('Number of full genes per insert')
    plt.savefig('full_genes_per_fragment.png', bbox_inches='tight')


def plot_partial_genes_per_fragment(insert_cov):
    """
    Plot the partial genes per insert
    :param insert_cov: DataFrame containing the number of partial genes per insert
    """
    if insert_cov.empty:
        plt.savefig('partial_genes_per_fragment.png') # empty plot
        write_empty_file('partial_genes_per_fragment.csv')
        return
    insert_cov.to_csv('partial_genes_per_fragment.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(insert_cov,
                 x='count',
                 hue='sample',
                 stat='frequency',
                 common_norm=True,
                 ax=ax)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_title('Number of partial or full genes per insert')
    plt.savefig('partial_genes_per_fragment.png', bbox_inches='tight')


def plot_barcode_length_boxplot(barcode_df):
    """
    Plot a boxplot of barcode lengths
    :param barcode_df: DataFrame containing barcode lengths
    """
    df_plot = barcode_df
    fig, ax = plt.subplots(figsize=(10, 6))
    if df_plot.empty:
        plt.savefig('barcode_length_distribution.png')  # empty plot
        write_empty_file('barcode_lengths.csv')
        return
    df_plot.to_csv('barcode_lengths.csv', index=False)

    sns.boxplot(data=df_plot, x="sample", y="barcode_len", ax=ax)
    plt.xticks(rotation=90)
    ax.set_title('Distribution of barcode lengths')
    plt.savefig('barcode_length_distribution.png', bbox_inches='tight')


def plot_proportions_of_barcodes(barcode_df):
    """
    Plot the proportion of true/empty/other detected barcodes
    :param barcode_df: DataFrame containing barcode information
    """
    if barcode_df.empty:
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
    """
    Plot the copy number of each unique barcode
    :param barcode_df: DataFrame containing barcode information
    """
    if barcode_df.empty:
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



def plot_insert_length_histogram(inserts):
    """
    Plot the distribution of insert lengths
    :param inserts: DataFrame containing insert information
    """
    df_plot = inserts
    if df_plot.empty:
        plt.savefig('insert_length_distribution.png') # empty plot
        write_empty_file('insert_length_distribution.csv')
        return
    df_plot.to_csv('insert_length_distribution.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.violinplot(df_plot,
                   x='insert_len',
                   y='sample',
                   ax=ax)
    ax.set_title('Distribution of insert lengths')
    plt.savefig('insert_length_distribution.png', bbox_inches='tight')

def concatenate_files(file_map, summary_type, output_file, save_file=True):
    """
    Given a file map, concatenate the files and return the concatenated DataFrame
    :param file_map: str representing the file map path
    :param summary_type: str representing the type of summary
    :param output_file: str representing the output file path
    :param save_file: bool representing whether to save the file
    :return: The concatenated DataFrame
    """
    with open(file_map, 'r') as f:
        # Create a dictionary of sample names and file paths
        sample_file_dict = {line.split('\t')[0]: line.split('\t')[1].rstrip('\n') for line in f}
    df_to_concatenate = []
    df = None
    for key,val in sample_file_dict.items():
        try:
            # Read in the data and, depending on the summary type, manipulate as needed
            if summary_type == 'barcode':
                print(f'File {val} sample {key}')
                df = pd.read_table(val, names=['read', 'barcode_seq', 'barcode_len'], usecols=[0, 1, 3], quoting=csv.QUOTE_NONE, engine='c', sep="\t")

            if summary_type == 'insert':
                print(f'Processing {val} for sample {key}')
                df = pd.read_table(val, names=['read', 'insert_seq', 'insert_len'], usecols=[0, 1, 3], engine='c', quoting=csv.QUOTE_NONE, sep="\t")

            if summary_type == 'insert_coverage':
                cov = pd.read_table(val, header=None, engine='c', sep="\t")
                df = cov.loc[:, [3, 6, 7, 8, 9]].rename(
                    columns={3: 'read', 6: 'count', 7: 'bases', 8: 'read_length', 9: 'percent_cov'})

            if summary_type == 'seq_stat':
                df = pd.read_table(val)

            if summary_type == 'vec_map':
                df = pd.read_table(val, names = ['value', 'key'], usecols=[0,2], engine='c', sep="\t")

            if summary_type == 'barcode_counts':
                df = pd.read_table(val, names=['barcode_seq', 'barcode_count'], engine='c', sep="\t")

            df['sample'] = key
            df_to_concatenate.append(df)
        except pd.errors.EmptyDataError:
            print(f'No seq stat data for {key} found in {val}')
            continue

    if not df_to_concatenate:
        write_empty_file(output_file)
        return pd.DataFrame() # return empty DataFrame

    concatenated_df = pd.concat(df_to_concatenate)
    if save_file:
        concatenated_df.to_csv(output_file, index=False)
    return concatenated_df


def seq_summary(barcode_data, insert_data, seq_stat, vec_map_stats, samp_info=None):
    """
    Generate a sequence summary table
    :param barcode_data: DataFrame containing barcode data
    :param insert_data: DataFrame containing insert data
    :param seq_stat: DataFrame containing sequence statistics
    :param vec_map_stats: DataFrame containing vector mapping statistics
    :param samp_info: DataFrame containing sample information
    """
    both = barcode_data.merge(insert_data, on=['sample', 'read']) \
        .groupby('sample', as_index=False) \
        .size().rename(columns={'size': 'both'})

    map_stats = vec_map_stats \
        .pivot(columns='key', values='value', index='sample').reset_index() \
        [['sample', 'primary mapped']]

    seq_stat["name"] = seq_stat["file"].apply(lambda x: str(Path(x).with_suffix('')))
    seq_stat["name"] = seq_stat["name"] \
        .apply(lambda x: x if x in ['inserts', 'barcodes'] else "raw_reads")
    num_seqs = seq_stat[['sample', 'name', 'num_seqs']] \
        .pivot(columns="name", values="num_seqs", index="sample") \
        .merge(both, on='sample', how='left') \
        .merge(map_stats) \
        [['sample', 'raw_reads', 'primary mapped', 'barcodes', 'inserts', 'both']] \
        .assign(pct_mapped=lambda x: round(100 * (x['primary mapped'].astype(int) / x.raw_reads), 2)) \
        .assign(bc_pct=lambda x: round(100 * (x.barcodes / x.raw_reads), 2)) \
        .assign(ins_pct=lambda x: round(100 * (x.inserts / x.raw_reads), 2)) \
        .assign(both_pct=lambda x: round(100 * (x.both / x.raw_reads), 2)) \
        .fillna(0) \
        .rename(columns={'raw_reads': 'Raw Reads',
                         'primary mapped': 'Reads mapped to plasmid',
                         'pct_mapped': '% mapped to plasmid',
                         'barcodes': 'Reads with barcodes', 'inserts': 'Reads with inserts',
                         'both': 'Reads with both', 'bc_pct': '% with barcode', 'ins_pct': '% with insert',
                         'both_pct': '% with both', 'Sample': 'sample'})

    if samp_info is not None:
        num_seqs = samp_info.merge(num_seqs, on='sample')

    return num_seqs


def process(sample_file_map, summary_type, **kwargs):
    """
    Process the data and generate the appropriate plots
    :param sample_file_map: str representing the sample file map
    :param summary_type: str representing the type of summary
    :param kwargs: dict containing additional arguments
    """

    match summary_type:
        case 'barcode':
            output_file_name = 'concatenated_barcodes.csv'
            concatenated_df = concatenate_files(sample_file_map, summary_type, output_file_name)
            plot_barcode_length_boxplot(concatenated_df)
            plot_proportions_of_barcodes(concatenated_df)
            plot_copy_number(concatenated_df)
        case 'insert':
            output_file_name = 'concatenated_inserts.csv'
            concatenated_df = concatenate_files(sample_file_map, summary_type, output_file_name)
            plot_insert_length_histogram(concatenated_df)
        case 'insert_coverage':
            output_file_name = 'concatenated_insert_coverage.csv'
            concatenated_df = concatenate_files(sample_file_map, summary_type, output_file_name)
            plot_full_genes_per_fragment(concatenated_df)
            plot_partial_genes_per_fragment(concatenated_df)
        case 'seq_stat':
            barcode_map = kwargs.get('barcode')
            insert_map = kwargs.get('insert')
            vector_map = kwargs.get('vector')
            seq_stats_df = concatenate_files(sample_file_map, summary_type, 'concatenated_seq_stats.csv')
            barcode_df = concatenate_files(barcode_map, 'barcode',
                                           None, False)
            insert_df = concatenate_files(insert_map, 'insert', None, False)
            vector_output_file_name = 'concatenated_vector_map_stats.csv'
            vector_df = concatenate_files(vector_map, 'vec_map', vector_output_file_name)
            # Combine all if possible
            if barcode_df.empty  or insert_df.empty  or seq_stats_df.empty or vector_df.empty:
                write_empty_file('seq_summary.csv')
                print('No seq data to summarize')
            else:
                num_seqs_df = seq_summary(barcode_df, insert_df, seq_stats_df, vector_df)
                num_seqs_df.to_csv('seq_summary.csv', index=False)

        case 'barcode_counts':
            output_file_name = 'concatenated_barcode_counts.csv'
            concatenate_files(sample_file_map, summary_type, output_file_name)

if __name__ == '__main__':
    if len(sys.argv) == 3:
        process(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 6:
        process(sys.argv[1], sys.argv[2], barcode=sys.argv[3], vector=sys.argv[4], insert=sys.argv[5])

