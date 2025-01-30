#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.io as pio
import os
import sys
from report_utils_template import load_report_data, seq_summary


sns.set_theme()
sns.set_style('darkgrid')
pioneer_colors = ['#FF8633', '#423759', '#314942', '#FFA632', '#F7F3ED']
sns.set_palette(sns.color_palette(pioneer_colors))

pio.templates['pioneer'] = pio.templates["seaborn"]
pio.templates['pioneer'].layout.colorway = pioneer_colors
pio.templates.default = 'pioneer'


def load_data(samples_path, result_dir):
    # sample data
    samples = pd.read_csv(samples_path)
    samps = {x: os.path.join(result_dir, x) for x in samples.id.to_list()}

    # load results
    data = load_report_data(samps)

    return data


def plot_barcode_length_histogram(data):
    df_plot = data.barcodes
    fig, ax = plt.subplots(figsize=(4, 4))
    if df_plot is None:
        plt.savefig('barcode_length_distribution.png') # empty plot
    sns.violinplot(df_plot,
                   x='barcode_len',
                   y='sample',
                   ax=ax)
    ax.set_title('Distribution of barcode lengths')
    plt.savefig('barcode_length_distribution.png')


def plot_proportions_of_barcodes(data):
    # Categorize barcodes by length and plot fraction composition for each library
    # Use 1bp tolerance on either end
    barcode_length = 44  # expected bp length of barcodes
    empty_cutoff = 5  # cutoff for labeling a vector as empty
    data.barcodes['barcode_exists'] = 'other'
    data.barcodes.loc[
        data.barcodes.barcode_len.between(barcode_length - 1, barcode_length + 1), 'barcode_exists'] = 'true'
    data.barcodes.loc[data.barcodes.barcode_len <= empty_cutoff, 'barcode_exists'] = 'false'

    df_plot = data.barcodes.groupby('sample')['barcode_exists'].value_counts(normalize=True).to_frame()
    fig, ax = plt.subplots(figsize=(6, 4))

    sns.barplot(df_plot,
                y='sample',
                x='proportion',
                hue='barcode_exists',
                ax=ax)
    ax.set_title('Proportion of true/empty/other detected barcodes')
    plt.savefig( 'barcode_proportions.png')


def plot_copy_number(data):
    df_plot = data.barcodes.groupby('sample')['barcode_seq'].value_counts().to_frame().rename(
        columns={'count': 'counts_per_barcode'})
    fig, ax = plt.subplots(figsize=(6, 4))
    if df_plot is None:
        plt.savefig('barcode_copy_number.png') # empty plot
    sns.histplot(df_plot,
                 x='counts_per_barcode',
                 hue='sample',
                 ax=ax)
    ax.set_title('Copy number of each unique barcode')
    plt.savefig( 'barcode_copy_number.png')


def plot_insert_length_histogram(data):
    df_plot = data.inserts
    if df_plot is None:
        plt.savefig('insert_length_distribution.png') # empty plot
    fig, ax = plt.subplots(figsize=(6, 4))

    sns.violinplot(df_plot,
                   x='insert_len',
                   y='sample',
                   ax=ax)
    ax.set_title('Distribution of insert lengths')
    plt.savefig('insert_length_distribution.png')


def plot_full_genes_per_fragment(data):
    if data.insert_cov_full is None:
        plt.savefig('full_genes_per_fragment.png') # empty plot
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.histplot(data.insert_cov_full,
                 x='count',
                 hue='sample',
                 stat='frequency',
                 common_norm=True,
                 ax=ax)
    ax.set_title('Number of full genes per insert')
    plt.savefig('full_genes_per_fragment.png')


def plot_partial_genes_per_fragment(data):
    if data.insert_cov is None:
        plt.savefig('partial_genes_per_fragment.png') # empty plot
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.histplot(data.insert_cov,
                 x='count',
                 hue='sample',
                 stat='frequency',
                 common_norm=True,
                 ax=ax)
    ax.set_title('Number of partial or full genes per insert')
    plt.savefig('partial_genes_per_fragment.png')


def visualize_results(samples_path, result_dir):
    # Deal with symlink
    true_path = os.readlink(samples_path)

    data = load_data(true_path, result_dir)

    # Export read statistics for raw reads/barcodes/inserts
    data.seq_stat.to_csv('raw_seq_stats.csv', index=False)


    # sequence summary
    num_seqs = seq_summary(data.barcodes, data.inserts, data.seq_stat, data.vec_map_stats)
    num_seqs.to_csv('seq_summary.csv', index=False)

    # plot barcode length histogram
    plot_barcode_length_histogram(data)

    # plot proportions of barcodes
    plot_proportions_of_barcodes(data)

    # plot copy number
    plot_copy_number(data)

    # plot insert length histogram
    plot_insert_length_histogram(data)

    # plot full genes per fragment
    plot_full_genes_per_fragment(data)

    # plot partial genes per fragment
    plot_partial_genes_per_fragment(data)



if __name__ == '__main__':
    visualize_results(sys.argv[1], sys.argv[2])



