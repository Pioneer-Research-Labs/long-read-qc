#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

def extract_coords(df_list, up_name, dn_name, filter_empty = True):
    up = df_list[up_name]
    dn = df_list[dn_name]

    # choose right column base on orientation
    up['start'] = np.where(up[[4]] == '+', up[[3]], up[[2]]) 
    up = up.rename(columns = {0: 'read'}).loc[:, ['read', 'start']]
    
    dn['end'] = np.where(dn[[4]] == '+', dn[[2]], dn[[3]])
    dn = dn.rename(columns = {0: 'read', 4: 'strand'}).loc[:, ['read', 'end', 'strand']]

    bed = up.merge(dn, on = 'read')
    bed['name'] = bed['read']
    bed['score'] = 0
    out = bed.loc[:, ['read', 'start', 'end', 'name', 'score', 'strand']]
    if filter_empty:
        nempty = len(out[out.start == out.end])
        print(f'Removing {nempty} empty reads for {up_name} and {dn_name}.')
        out = out[out.start != out.end]
    return(out)


paf_file = sys.argv[1]

maps = pd.read_table(paf_file, header  = None)

h = [g for _, g in maps.groupby(5)]
map_list = {x[5].unique().tolist()[0]:x for x in h}


ins_bed = extract_coords(map_list, 'INSERTUP', 'INSERTDN')
bc_bed = extract_coords(map_list, 'BARCODEUP', 'BARCODEDN')


bc_bed.to_csv("barcode_coords.bed", sep = "\t", index = False, header = False)
ins_bed.to_csv("insert_coords.bed", sep = "\t", index = False, header = False)
