#!/usr/bin/env python

import pandas as pd
import numpy as np
import json
import sys


def filter_bed(path, out_bed, out_stats):
    bed = pd.read_table(path, header=None)
    stats = {}
    stats['total_hits'] = bed.shape[0]

    hit_counts = bed[0].value_counts() 
    dup_list = hit_counts[hit_counts > 1].index.tolist()
    stats['n_duplicate_hits'] = len(dup_list)

    # keep only reads with one hit
    no_dups = bed[~bed[0].isin(dup_list)]

    # How many 0 size barcodes/inserts
    zero_size = no_dups[1] >= no_dups[2]
    stats['n_zero_size'] = zero_size.values.sum().item()

    out = no_dups[~zero_size]

    stats['n_good'] = out.shape[0]

    with open(out_stats, 'w') as fp:
        json.dump(stats, fp)

    out.to_csv(out_bed, sep = "\t", index = False, header = False)



filter_bed(sys.argv[1], sys.argv[2], sys.argv[3])