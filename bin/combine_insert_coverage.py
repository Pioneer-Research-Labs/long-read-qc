#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def write_empty_file(filename):
    with open(filename, "w") as my_empty_csv:
        pass


def plot_full_genes_per_fragment(insert_cov_full):
    if insert_cov_full is None:
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
    if insert_cov is None:
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


def concatenate_insert_coverage_files(file_map):
    with open(file_map, 'r') as f:
        sample_file_dict = {line.split()[0]: line.split()[1] for line in f}
    df_to_concatenate = []
    for key,val in sample_file_dict.items():
        insert_coverage_df = pd.read_table(val, header=None)
        read_cov = insert_coverage_df.loc[:, [3, 6, 7, 8, 9]].rename(
            columns={3: 'read', 6: 'count', 7: 'bases', 8: 'read_length', 9: 'percent_cov'})
        read_cov['sample'] = key
        df_to_concatenate.append(read_cov)

    all_insert_coverage_df = pd.concat(df_to_concatenate)
    all_insert_coverage_df.to_csv('concatenated_insert_coverage.csv', index=False)
    return all_insert_coverage_df

if __name__ == '__main__':
    # Load insert coverage from all samples
    coverage_df = concatenate_insert_coverage_files(sys.argv[1])
    # Plot full gene coverage
    plot_full_genes_per_fragment(coverage_df)
    # Plot partial gene coverage
    plot_partial_genes_per_fragment(coverage_df)
