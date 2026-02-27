#!/usr/bin/env python

import pandas as pd
import sys
import threading
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import os
import plotly.io as pio


# Color settings for seaborn and plotly
sns.set_theme(font_scale=.8)
sns.set_style('darkgrid')
pioneer_colors = ['#FF8633', '#423759', '#314942', '#FFA632', '#F7F3ED']
sns.set_palette(sns.color_palette(pioneer_colors))

pio.templates['pioneer'] = pio.templates["seaborn"]
pio.templates['pioneer'].layout.colorway = pioneer_colors
pio.templates.default = 'pioneer'


def process_genome_tags(cutadapt_file, tesseract_oligos, sample_name):
    """
    Process genome tags from a CSV file and save the results to a new CSV file.

    Parameters:
    input_file (str): Path to the input CSV file containing genome tags.
    """
    # Load the small CSV file into memory
    tess_df = pd.read_csv(tesseract_oligos, usecols=['8bp Barcode','Base Strain'])
    # Upper case the '8bp Barcode' column to ensure consistency
    tess_df['8bp Barcode'] = tess_df['8bp Barcode'].str.upper()
    tess_df['Base Strain'] = tess_df['Base Strain'].str.replace(" ", "_")
    df = tess_df[tess_df['Base Strain'].notna()]

    # Prepare an empty list to collect merged chunks
    merged_chunks = []

    # Read the large TSV file in chunks
    chunk_iter = pd.read_csv(cutadapt_file,names=['seqid','8bp Barcode'], header=None, sep='\t', chunksize=100000)

    for chunk in chunk_iter:
        merged = chunk.merge(df, on='8bp Barcode')  # replace 'common_column' with your join key
        merged_chunks.append(merged)

    # Concatenate all merged chunks into a single DataFrame
    genome_tags = pd.concat(merged_chunks, ignore_index=True)
    # Save the merged DataFrame to a new CSV file with columns 'seqid', 'Base Strain', and 'Genome Path'
    genome_tags.to_csv(sample_name + "_assigned_genome_tags.csv", index=False)

    file_outputs = []
    unique_genomes = genome_tags['Base Strain'].unique()
    for genome in unique_genomes:
        genome_df = genome_tags[genome_tags['Base Strain'] == genome]
        # Remove " rc" from the 'seqid' column and save to a new CSV file with only the 'seqid' column, essentially
        # creating a file of sequences ids that will later be used to extract sequences from the fastq file using seqkit.
        # We're removing " rc" because some of the sequence ids in the cutadapt output have " rc" appended to them (to indicate reverse complement)
        # and we want to ensure that the sequence ids in the output files match the sequence ids in the fastq file.
        genome_df_clean = genome_df['seqid'].str.replace(" rc", "")
        output_filename = f"{genome}.csv"
        genome_df_clean.to_csv(output_filename,columns=['seqid'], index=False,header=False)
        file_outputs.append(output_filename)
    return file_outputs


def run_command(cmd):
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        print(f"Output: {out.decode()}")
        if err:
                print(f"Error: {err.decode()}")


def run_seqkit(output_files, fastq_file, sample_name, construct, outdir, tesseract_oligo_grid):
    """
    Run seqkit to extract sequences from a fastq.gz file using the output CSV files.

    Parameters:

    output_files (list): List of output CSV files containing genome sequences.
    """
    commands = []
    for file in output_files:
        # grep -f id.txt seqs.fq.gz -o result.fq.gz
        cmd = f"seqkit grep -f {file} {fastq_file} -o {file.replace('.csv', '.fq.gz')}"
        commands.append(cmd)

    threads = []
    for cmd in commands:
        t = threading.Thread(target=run_command, args=(cmd,))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()

    # Generate a dictionary from the tesseract oligo grid to map the genome names to the fully qualified genome names
    # found on s3 in pioneer-genomes
    genome_dict = convert_file_to_dict(tesseract_oligo_grid)
    with open(f"{sample_name}_sample_sheet.csv", 'w') as f:
        fastq_files = []
        for file in output_files:
            fastq = file.replace('.csv', '.fq.gz')
            fastq_files.append(fastq)
            genome = file.split('.')[0]
            s3_location = f"{outdir}/{sample_name}/genome_fastqs/{fastq}"
            genome_path = genome_dict[genome]  # Get the fully qualified genome name from the dictionary, or use the original genome name if not found
            f.write(f"{sample_name}_{genome},{genome_path},{construct},{s3_location}\n")

    return fastq_files

def convert_file_to_dict(file_name, key_col=2, value_col=3):
    """
    Convert a CSV file to a dictionary given the key and value CSV columns.

    Parameters:
    file_name (str): Path to the input CSV file.

    Returns:
    dict: A dictionary with keys from the first column and values from the second column.
    """
    df = pd.read_csv(file_name)
    return dict(zip(df.iloc[:, key_col], df.iloc[:, value_col]))


def generate_seq_summaries(fastqs):
    """
    Generate sequence summaries for the extracted sequences using seqkit stats.
    """
    # Take all the fastqs in the current directory that end with .fq.gz and run seqkit stats on them,
    # saving the output to a new file called "genome_name_seq_summary.txt"
    for fastq in fastqs:
        genome_name = fastq.split('.')[0]
        cmd = f"seqkit stats {fastq} --tabular > {genome_name}_seq_summary.txt"
        run_command(cmd)

    # Combine all the individual genome sequence summary files into a single file called "combined_seq_summary.txt"
    with open("combined_seq_summary.txt", 'w') as outfile:
        for fastq in fastqs:
            genome_name = fastq.split('.')[0]
            summary_file = f"{genome_name}_seq_summary.txt"
            with open(summary_file) as infile:
                #Write the second line of the summary file (which contains the actual stats)
                # to the combined summary file, along with the genome name for reference
                lines = infile.readlines()
                if len(lines) > 1:
                    stats = lines[1]
                    outfile.write(stats) # Add a newline between summaries
    # Add headers to the combined summary file
    combined_df = pd.read_csv("combined_seq_summary.txt",
                              sep='\t',
                              names=["file", "format", "type" , "num_seqs", "sum_len", "min_len", "avg_len","max_len"] )
    combined_df.to_csv("combined_seq_summary.tsv", index=False, sep='\t')
    return combined_df

def generate_barplot(combined_seq_summary_df):
    """
    Generate a barplot of sequence count per genome.
    """
    # Create a new row in the dataframe that subtracts all rows from second column in
    # the last row (which is the original fastq file) to get the number of sequences that were not assigned to any genome, and label this row "unassigned"
    unassigned_num_seqs = combined_seq_summary_df['num_seqs'].iloc[-1] - combined_seq_summary_df['num_seqs'].iloc[:-1].sum()
    unassigned_row = pd.DataFrame({'file': ['unassigned reads'], 'format': [''], 'type': [''], 'num_seqs': [unassigned_num_seqs], 'sum_len': [''], 'min_len': [''], 'avg_len': [''], 'max_len': ['']})
    combined_seq_summary_df = pd.concat([combined_seq_summary_df, unassigned_row], ignore_index=True)
    # Order the combined summary dataframe by the number of sequences in descending order
    combined_seq_summary_df = combined_seq_summary_df.sort_values(by='num_seqs', ascending=False)
    # Strip the fq.gz from the file names in the 'file' column for better visualization
    combined_seq_summary_df['file'] = combined_seq_summary_df['file'].str.replace('.fq.gz', '', regex=False)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='file', y='num_seqs', data=combined_seq_summary_df)
    plt.xticks(rotation=90)
    plt.xlabel('Genome')
    plt.ylabel('Number of Sequences')
    plt.title('Number of Sequences Found with 8bp Barcodes by Genome')
    plt.tight_layout()
    plt.savefig('num_sequences_per_genome.png')
    plt.close()



if __name__ == '__main__':
    genome_tags_tsv = sys.argv[1]
    tesseract_oligos_csv = sys.argv[2]
    sample_name = sys.argv[3]
    fastq_file = sys.argv[4]
    construct = sys.argv[5]
    output_dir: str = sys.argv[6]
    outputs = process_genome_tags(genome_tags_tsv, tesseract_oligos_csv, sample_name)
    seq_fq_files = run_seqkit(outputs, fastq_file, sample_name, construct, output_dir, tesseract_oligos_csv)
    new_name = f"original_{fastq_file}"
    os.rename(fastq_file, new_name)  # Rename the original fastq file to keep it for reference
    seq_fq_files.append(new_name)

    combined_summary = generate_seq_summaries(seq_fq_files)
    generate_barplot(combined_summary)

