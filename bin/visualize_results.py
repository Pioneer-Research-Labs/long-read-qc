#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.io as pio
import os
import sys
from report_utils_template import load_report_data, seq_summary


sns.set_theme()
sns.set(font_scale=.8)
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

def write_empty_file(filename):
    with open(filename, "w") as my_empty_csv:
            pass

def plot_barcode_length_histogram(data):
    df_plot = data.barcodes
    fig, ax = plt.subplots(figsize=(10, 6))
    if df_plot is None:
        plt.savefig('barcode_length_distribution.png') # empty plot
        write_empty_file('barcode_lengths.csv')
        return
    df_plot.to_csv('barcode_lengths.csv', index=False)
    
    sns.boxplot(data=df_plot, x="sample", y="barcode_len", ax=ax)
    ax.set_title('Distribution of barcode lengths')
    plt.savefig('barcode_length_distribution.png', bbox_inches='tight')


def plot_proportions_of_barcodes(data):
    if data.barcodes is None:
        plt.savefig('barcode_proportions.png')
        write_empty_file('barcode_proportions.csv')
        return
    # Categorize barcodes by length and plot fraction composition for each library
    # Use 1bp tolerance on either end
    barcode_length = 44  # expected bp length of barcodes
    empty_cutoff = 5  # cutoff for labeling a vector as empty
    data.barcodes['barcode_exists'] = 'other'
    data.barcodes.loc[
        data.barcodes.barcode_len.between(barcode_length - 1, barcode_length + 1), 'barcode_exists'] = 'true'
    data.barcodes.loc[data.barcodes.barcode_len <= empty_cutoff, 'barcode_exists'] = 'false'

    df_plot = data.barcodes.groupby('sample')['barcode_exists'].value_counts(normalize=True).to_frame()
    data.barcodes.to_csv('barcode_proportions.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.barplot(df_plot,
                y='sample',
                x='proportion',
                hue='barcode_exists',
                ax=ax)
    ax.set_title('Proportion of true/empty/other detected barcodes')
    plt.savefig( 'barcode_proportions.png', bbox_inches='tight')


def plot_copy_number(data):
    if data.barcodes is None:
        plt.savefig('barcode_copy_number.png') # empty plot
        write_empty_file('barcode_copy_number.csv')
        return
    df_plot = data.barcodes.groupby('sample')['barcode_seq'].value_counts().to_frame().rename(
        columns={'count': 'counts_per_barcode'})
    
    fig, ax = plt.subplots(figsize=(10, 6))

    df_plot.to_csv('barcode_copy_number.csv', index=False)
    sns.histplot(df_plot,
                 x='counts_per_barcode',
                 hue='sample',
                 ax=ax)
    ax.set_title('Copy number of each unique barcode')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig( 'barcode_copy_number.png', bbox_inches='tight')


def plot_insert_length_histogram(data):
    df_plot = data.inserts
    if df_plot is None:
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


def plot_full_genes_per_fragment(data):
    if data.insert_cov_full is None:
        plt.savefig('full_genes_per_fragment.png') # empty plot
        write_empty_file('full_genes_per_fragment.csv')
        return
    data.insert_cov_full.to_csv('full_genes_per_fragment.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(data.insert_cov_full,
                 x='count',
                 hue='sample',
                 stat='frequency',
                 common_norm=True,
                 ax=ax)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_title('Number of full genes per insert')
    plt.savefig('full_genes_per_fragment.png', bbox_inches='tight')


def plot_partial_genes_per_fragment(data):
    if data.insert_cov is None:
        plt.savefig('partial_genes_per_fragment.png') # empty plot
        write_empty_file('partial_genes_per_fragment.csv')
        return
    data.insert_cov.to_csv('partial_genes_per_fragment.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(data.insert_cov,
                 x='count',
                 hue='sample',
                 stat='frequency',
                 common_norm=True,
                 ax=ax)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_title('Number of partial or full genes per insert')
    plt.savefig('partial_genes_per_fragment.png', bbox_inches='tight')


def visualize_results(samples_path, result_dir):

    data = load_data(samples_path, result_dir)

    # Export read statistics for raw reads/barcodes/inserts
    if data.seq_stat is None:
        write_empty_file('raw_seq_stats.csv')
    else:
        data.seq_stat.to_csv('raw_seq_stats.csv', index=False)


    # sequence summary
    if data.barcodes is None or data.inserts is None or data.seq_stat is None or data.vec_map_stats is None:
        write_empty_file('seq_summary.csv')

    else:
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
    c = os.getcwd().split("work")[0]
    absolute_path = os.path.join(c,"results",sys.argv[1])
    print(absolute_path)
    visualize_results(absolute_path, sys.argv[2])



