import pandas as pd
import json
import os
from pathlib import Path



def list_files(samps, basename):
    files = {key:os.path.join(val, basename) for key,val in samps.items()}
    return {key:val for key,val in files.items() if os.path.exists(val)}

def load_barcode_data(path, samp):
    try:
        barcodes = pd.read_table(path, names = ['read', 'barcode_seq', 'barcode_len'], usecols = [0,1,3])
        barcodes['sample'] = samp
    except:
        barcodes = None
    return barcodes

def load_barcode_counts(path, samp):
    try:
        counts = pd.read_table(path, names = ['barcode_seq', 'barcode_count']) \
            .assign(sample = samp)
    except:
        counts = None
    return counts

def load_insert_data(path, samp):
    try:
        inserts = pd.read_table(path, names = ['read', 'insert_seq', 'insert_len'],  usecols = [0,1,3])
        inserts['sample'] = samp
    except:
        inserts = None
    return inserts

def load_genome_cov(path, samp):
    genome_cov = pd.read_table(path, names = ['chr', 'pos', 'cov'])
    genome_cov['sample'] = samp
    return genome_cov

def load_gene_cov(path, samp):
    gene_cov_full = pd.read_table(path, header = None)
    gene_cov = gene_cov_full.loc[:, [9, 10, 11, 12]].rename(
        columns={9:'count', 10:'bases', 11:'gene_length', 12:'percent_cov'})
    gene_cov['sample'] = samp
    return gene_cov

def load_insert_cov(path, samp):
    try:
        read_cov_full = pd.read_table(path, header = None)
        read_cov = read_cov_full.loc[:, [3, 6, 7, 8, 9]].rename(
            columns={3: 'read', 6:'count', 7:'bases', 8:'read_length', 9:'percent_cov'})
        read_cov['sample'] = samp
    except:
        read_cov = None
    return read_cov
    
def load_seq_stats(path, samp):
    df = pd.read_table(path)
    df['sample'] = samp
    return(df)

def load_cov_stats(path, samp):
    df = pd.read_table(path)
    df['sample'] = samp
    return(df)

def read_div(n, d):
    return n / d if d else 0


def load_cutadapt_report(samps, basename):

    cut_reports = {}
    for key,val in list_files(samps, basename).items():
        with open(val, 'r') as jf:
            cut_reports[key] = (json.load(jf))

    report_list = []
    for key,val in cut_reports.items():
        in_reads = val['read_counts']['input']
        out_reads = val['read_counts']['output']
        discard = in_reads - out_reads
        frac = round(100*read_div(out_reads, in_reads), 2)
        frac_discard = round(100*read_div(discard, in_reads), 2)
        out = {'Sample': [key], 'Input Reads': [in_reads], 'Num detected': [out_reads], 
            'Num discarded': [discard], '% detected': [frac], '% discarded': [frac_discard]}
        report_list.append(pd.DataFrame.from_dict(out))

    return pd.concat(report_list)

def load_intersects(path, samp):
        cols = ['ins_chr', 'ins_start', 'ins_end', 'ins_id', 'ins_score', 'ins_strand',
                'gene_chr', 'gene_start', 'gene_end', 'gene_name', 'gene_score', 'gene_strand', 'overlap']
        insert_ovlps = pd.read_table(path, names = cols)[['ins_id', 'ins_start', 'ins_end', 'gene_name', 
                                                                'gene_start', 'gene_end', 'gene_strand', 'overlap']]
        insert_ovlps.insert(0, 'sample', samp)

        out = insert_ovlps \
                .assign(gene_size = lambda x: x.gene_end - x.gene_start) \
                .assign(ins_size = lambda x: x.ins_end - x.ins_start) \
                .assign(pct_gene_ovlp = lambda x: 100 * (x.overlap / x.gene_size)) \
                .assign(pct_ins_ovlp = lambda x: 100 * (x.overlap / x.ins_size)) \
                .sort_values('ins_start')                               
        return out

def seq_summary(barcode_data, insert_data, seq_stat, samp_info = None):
    both = barcode_data.merge(insert_data, on = ['sample', 'read']) \
        .groupby('sample', as_index = False) \
        .size().rename(columns = {'size': 'both'})

    seq_stat["name"] = seq_stat["file"].apply(lambda x: str(Path(x).with_suffix('')))
    seq_stat["name"] = seq_stat["name"] \
        .apply(lambda x: x if x in ['reads_rotated', 'inserts', 'barcodes'] else "raw_reads")
    num_seqs = seq_stat[['sample', 'name', 'num_seqs']] \
        .pivot(columns = "name", values = "num_seqs", index = "sample") \
        .merge(both, on = 'sample', how = 'left') \
        [['sample', 'raw_reads', 'reads_rotated', 'barcodes', 'inserts', 'both']] \
        .assign(pl_pct = lambda x: round(100*(x.reads_rotated / x.raw_reads), 2)) \
        .assign(bc_pct = lambda x: round(100*(x.barcodes / x.raw_reads), 2)) \
        .assign(ins_pct = lambda x: round(100*(x.inserts / x.raw_reads), 2)) \
        .assign(both_pct = lambda x: round(100*(x.both / x.raw_reads), 2)) \
        .fillna(0) \
        .rename(columns = {'raw_reads': 'Raw Reads', 'reads_rotated': "Reads with plasmid sequence", 
                        'barcodes': 'Reads with barcodes', 'inserts': 'Reads with inserts',
                        'both': 'Reads with both', 'bc_pct': '% with barcode', 'ins_pct': '% with insert',
                        'pl_pct': '% with plasmid', 'both_pct': '% with both', 'Sample': 'sample'})

    if samp_info is not None:
        num_seqs = samp_info.merge(num_seqs, on = 'sample')
        
    return(num_seqs)