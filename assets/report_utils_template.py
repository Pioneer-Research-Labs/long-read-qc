import pandas as pd
import json
import os
from pathlib import Path



# for organizing our data a bit better
class PipelineData(dict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as k:
            raise AttributeError(k)

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError as k:
            raise AttributeError(k)

    def __repr__(self):
        return '<PipelineData' + dict.__repr__(self) + '>'


def load_report_data(samps):
    out = PipelineData({x['name']: load_sample_data(samps, x['load_fun'], x['path']) for x in to_load})
    return out

def list_files(samps, basename):
    files = {key:os.path.join(val, basename) for key,val in samps.items()}
    return {key:val for key,val in files.items() if os.path.exists(val)}

def load_sample_data(samps, load_fun, path):
    try:
        out = pd.concat([load_fun(val, key) for key, val in list_files(samps, path).items()])
    except:
        out = None
    return(out)


def load_csv(path, samp):
    try:
        x = pd.read_csv(path)
        x['sample'] = samp
    except:
        x = None
    return x 

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
    try:
        genome_cov = pd.read_table(path, names = ['chr', 'pos', 'cov'])
        genome_cov['sample'] = samp
    except:
        genome_cov = None
    return genome_cov


def load_map_stats(path, samp):
    try:
        x = pd.read_table(path, names = ['value', 'key'], usecols=[0,2])
        x['sample'] = samp
    except:
        x = None
    return(x)

def load_gene_cov(path, samp):
    try:
        gene_cov_full = pd.read_table(path, header = None)
        gene_cov = gene_cov_full.loc[:, [9, 10, 11, 12]].rename(
            columns={9:'count', 10:'bases', 11:'gene_length', 12:'percent_cov'})
        gene_cov['sample'] = samp
    except:
        gene_cov = None
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

def seq_summary(barcode_data, insert_data, seq_stat, vec_map_stats, samp_info = None):
    both = barcode_data.merge(insert_data, on = ['sample', 'read']) \
        .groupby('sample', as_index = False) \
        .size().rename(columns = {'size': 'both'})
    
    map_stats = vec_map_stats \
        .pivot(columns='key', values='value', index='sample').reset_index() \
        [['sample', 'primary mapped']] 

    seq_stat["name"] = seq_stat["file"].apply(lambda x: str(Path(x).with_suffix('')))
    seq_stat["name"] = seq_stat["name"] \
        .apply(lambda x: x if x in ['inserts', 'barcodes'] else "raw_reads")
    num_seqs = seq_stat[['sample', 'name', 'num_seqs']] \
        .pivot(columns = "name", values = "num_seqs", index = "sample") \
        .merge(both, on = 'sample', how = 'left') \
        .merge(map_stats) \
        [['sample', 'raw_reads', 'primary mapped', 'barcodes', 'inserts', 'both']] \
        .assign(pct_mapped = lambda x: round(100*(x['primary mapped'].astype(int) / x.raw_reads), 2)) \
        .assign(bc_pct = lambda x: round(100*(x.barcodes / x.raw_reads), 2)) \
        .assign(ins_pct = lambda x: round(100*(x.inserts / x.raw_reads), 2)) \
        .assign(both_pct = lambda x: round(100*(x.both / x.raw_reads), 2)) \
        .fillna(0) \
        .rename(columns = {'raw_reads': 'Raw Reads',  
                           'primary mapped': 'Reads mapped to plasmid',
                           'pct_mapped': '% mapped to plasmid',
                            'barcodes': 'Reads with barcodes', 'inserts': 'Reads with inserts',
                            'both': 'Reads with both', 'bc_pct': '% with barcode', 'ins_pct': '% with insert',
                            'both_pct': '% with both', 'Sample': 'sample'})

    if samp_info is not None:
        num_seqs = samp_info.merge(num_seqs, on = 'sample')
        
    return(num_seqs)


to_load = [
    {'name': 'barcodes', 'load_fun': load_barcode_data, 'path': 'barcodes.tsv'},
    {'name': 'barcode_counts', 'load_fun': load_barcode_counts, 'path': 'barcode_counts.tsv'},
    {'name': 'inserts', 'load_fun': load_insert_data, 'path': 'inserts.tsv'},
    {'name': 'genome_cov', 'load_fun': load_genome_cov, 'path': 'genome_coverage.tsv'},
    {'name': 'gene_cov', 'load_fun': load_gene_cov, 'path': 'gene_coverage.bed'},
    {'name': 'insert_cov', 'load_fun': load_insert_cov, 'path': 'insert_coverage.bed'},
    {'name': 'insert_cov_full', 'load_fun': load_insert_cov, 'path': 'insert_coverage_full.bed'},
    {'name': 'intersect', 'load_fun': load_intersects, 'path': 'insert_intersect.out'},
    {'name': 'seq_stat', 'load_fun': load_seq_stats, 'path': 'seq_stats.tsv'},
    {'name': 'vec_map_stats', 'load_fun': load_map_stats, 'path': 'mapped_vector_stats.tsv'},
    {'name': 'ins_map_stats', 'load_fun': load_map_stats, 'path': 'mapped_insert_stats.tsv'},
    {'name': 'cov_stat', 'load_fun': load_cov_stats, 'path': 'genome_cov_stats.tsv'},
    {'name': 'matches', 'load_fun': load_csv, 'path': 'insert_matches.csv'},
    {'name': 'tax', 'load_fun': load_csv, 'path': 'insert_taxonomy.csv'}
]


# taxonomy

def rank_summary(tax, rank):
    
    all_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    i = [all_ranks.index(r) + 1 for r in all_ranks if r  == rank]
    ranks = all_ranks[0:i[0]]

    x = tax.query(f'rank == "{rank}"')[['fraction', 'lineage', 'sample']]
    x[ranks] = x.lineage.str.split(';', expand = True)
    x = x.drop(['lineage'], axis = 1)
    x[ranks] = x[ranks].apply(lambda x: x.str.replace("(^[a-z]__)", "", regex = True)) 
    x['Percentage'] = round(x['fraction'] * 100, 2)

    return x

