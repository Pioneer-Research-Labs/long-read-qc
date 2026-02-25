#!/usr/bin/env python

import pandas as pd
import sys
import threading
import subprocess

def process_genome_tags(cutadapt_file, tesseract_oligos, sample_name):
    """
    Process genome tags from a CSV file and save the results to a new CSV file.

    Parameters:
    input_file (str): Path to the input CSV file containing genome tags.
    """
    # Load the small CSV file into memory
    tess_df = pd.read_csv(tesseract_oligos, usecols=['8bp Barcode','Base Strain', 'Genome Path'])
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

    output_files = []
    unique_genomes = genome_tags['Base Strain'].unique()
    for genome in unique_genomes:
        genome_df = genome_tags[genome_tags['Base Strain'] == genome]
        genome_df['seq_id'] = genome_df['seqid'].str.replace(" rc", "")
        output_filename = f"{genome}.csv"
        genome_df.to_csv(output_filename,columns=['seqid', 'Genome Path'], index=False,header=False)
        output_files.append(output_filename)
    return output_files


def run_command(cmd):
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        print(f"Output: {out.decode()}")
        if err:
                print(f"Error: {err.decode()}")


def run_seqkit(output_files, fastq_file, sample_name, construct, outdir):
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

    with open(f"{sample_name}_sample_sheet.csv", 'w') as f:
        for file in output_files:
            # The reference genome name (e.g. Bacillus_subtilis_PY79_GCF_023521615.1) is the second column in the output CSV file
            # We just need the first line of the output CSV file to get the reference genome name,
            # since all lines in the output CSV file have the same reference genome name
            with open(file, 'r') as g:
                first_line = g.readline()
                genome_path = first_line.split(',')[1]
                # Remove newline character from genome_path
                genome_path = genome_path.strip()
            fastq = file.replace('.csv', '.fq.gz')
            s3_location = f"{outdir}/{sample_name}/{fastq}"
            f.write(f"{sample_name},{genome_path},{construct},{s3_location}\n")

    return output_files

if __name__ == '__main__':
    genome_tags_tsv = sys.argv[1]
    tesseract_oligos_csv = sys.argv[2]
    sample_name = sys.argv[3]
    fastq_file = sys.argv[4]
    construct = sys.argv[5]
    output_dir: str = sys.argv[6]
    outputs = process_genome_tags(genome_tags_tsv, tesseract_oligos_csv, sample_name)
    output_files = run_seqkit(outputs, fastq_file, sample_name, construct, output_dir)

