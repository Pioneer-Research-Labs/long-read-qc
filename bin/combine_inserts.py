#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys


def write_empty_file(filename):
    with open(filename, "w") as my_empty_csv:
        pass


def plot_insert_length_histogram(inserts):
    df_plot = inserts
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


def concatenate_insert_files(file_map):
    with open(file_map, 'r') as f:
        sample_file_dict = {line.split()[0]: line.split()[1] for line in f}
    df_to_concatenate = []
    for key,val in sample_file_dict.items():
        inserts = pd.read_table(val, names=['read', 'insert_seq', 'insert_len'], usecols=[0, 1, 3])
        inserts['sample'] = key
        df_to_concatenate.append(inserts)

    insert_df = pd.concat(df_to_concatenate)
    insert_df.to_csv('concatenated_inserts.csv', index=False)
    return insert_df


if __name__ == '__main__':
    insert_map = sys.argv[1]
    insert_df = concatenate_insert_files(insert_map)
    plot_insert_length_histogram(insert_df)