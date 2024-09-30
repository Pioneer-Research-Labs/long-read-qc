import pandas as pd
import json
import os



def list_files(samps, basename):
    files = {key:os.path.join(val, basename) for key,val in samps.items()}
    return {key:val for key,val in files.items() if os.path.exists(val)}

def load_barcode_data(paths, samp):
    barcode_path, counts_path = paths
    barcodes = pd.read_table(barcode_path, names = ['read', 'barcode_seq', 'barcode_len'], usecols = [0,1,3])
    barcode_counts = pd.read_table(counts_path, names = ['barcode_seq', 'barcode_count'])
    barcodes = barcodes.merge(barcode_counts, on = 'barcode_seq')
    barcodes['sample'] = samp
    return barcodes

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

