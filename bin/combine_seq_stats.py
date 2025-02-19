#!/usr/bin/env python

import pandas as pd
from pathlib import Path
import sys
from combine_barcodes import concatenate_barcode_files
from combine_inserts import concatenate_insert_files


def write_empty_file(filename):
    with open(filename, "w") as my_empty_csv:
        pass


def concatenate_seq_stat_files(file_map):
    with open(file_map, 'r') as f:
        sample_file_dict = {line.split()[0]: line.split()[1] for line in f}
    df_to_concatenate = []
    for key,val in sample_file_dict.items():
        try:
            df = pd.read_table(val)
            df['sample'] = key
            df_to_concatenate.append(df)
        except pd.errors.EmptyDataError:
            print(f'No seq stat data for {key} found in {val}')
            continue

    if not df_to_concatenate:
        write_empty_file('concatenated_seq_stats.csv')
        return None
    seq_stats_df = pd.concat(df_to_concatenate)
    seq_stats_df.to_csv('concatenated_seq_stats.csv', index=False)
    return seq_stats_df


def seq_summary(barcode_data, insert_data, seq_stat, vec_map_stats, samp_info=None):
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


def concatenate_vec_map_stats_files(file_map):
    with open(file_map, 'r') as f:
        sample_file_dict = {line.split()[0]: line.split()[1] for line in f}
    df_to_concatenate = []
    for key, val in sample_file_dict.items():
        try:
            vec_stats = pd.read_table(val, names = ['value', 'key'], usecols=[0,2])
            vec_stats['sample'] = key
            df_to_concatenate.append(vec_stats)
        except pd.errors.EmptyDataError:
            print(f'No vector mapping stats for {key} found in {val}')
            continue

    if not df_to_concatenate:
        write_empty_file('concatenated_vector_map_stats.csv')
        return None
    vec_map_stats_df = pd.concat(df_to_concatenate)
    vec_map_stats_df.to_csv('concatenated_vector_map_stats.csv', index=False)
    return vec_map_stats_df

if __name__ == '__main__':
    # Load seq stats from all samples
    seq_stats = concatenate_seq_stat_files(sys.argv[1])
    # Load barcodes from all samples
    barcodes = concatenate_barcode_files(sys.argv[2])
    # Load vector mapping stats from all samples
    vec_map_stats = concatenate_vec_map_stats_files(sys.argv[3])
    # Load inserts from all samples
    inserts = concatenate_insert_files(sys.argv[4])

    # Combine all if possible
    if barcodes is None or inserts is None or seq_stats is None or vec_map_stats is None:
        num_seqs_df = seq_summary(barcodes, inserts, seq_stats, vec_map_stats)
        num_seqs_df.to_csv('seq_summary.csv', index=False)
    else:
        write_empty_file('seq_summary.csv')
        print('No seq data to summarize')