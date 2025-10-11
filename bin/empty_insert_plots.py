import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def empty_insert_plots(barcode_file, insert_file, site_file, title):
    barcodes = pd.read_csv(barcode_file, sep='\t', names=['id', 'seq', 'quality', 'length'])
    inserts = pd.read_csv(insert_file, sep='\t', names=['id', 'seq', 'quality', 'length'])
    sites = pd.read_csv(site_file, sep='\t', names=['id', 'seq', 'quality', 'length'])

    # Get site ids that have barcode but no inserts
    sites_with_barcodes = barcodes.merge(sites, on='id', copy=False)

    # Sites with barcodes and inserts
    sites_with_barcode_and_insert = sites_with_barcodes.merge(inserts, on='id', copy=False)

    # Sites with barcodes but no inserts
    sites_with_barcodes_no_inserts = sites_with_barcodes.merge(inserts, on='id', how='outer', indicator=True, ).query(
        '_merge=="left_only"').drop('_merge', axis=1)

    sns.histplot(data=sites_with_barcodes_no_inserts, x="length_x")
    plt.title(f"{title} - Histogram of sites with barcodes but no inserts")
    output_file = Path('histogram_of_sites_with_barcode_no_insert.png')
    plt.savefig(output_file)
    plt.close()

    # Write out file summarizing lengths for each dataframe
    with open("summary_of_counts.csv", "w") as summary_file:
        summary_file.write(
            f"Count of sequences with inserts,"
            f"Count of sequences with barcodes,"
            f"Count of sequences with site flanks,"
            f"Count of sequences with site flanks and barcodes,"
            f"Count of sequences with site flanks and barcodes but no inserts\n"
        )
        summary_file.write(
            f"{len(inserts)},"
            f"{len(barcodes)},"
            f"{len(sites)},"
            f"{len(sites_with_barcode_and_insert)},"
            f"{len(sites_with_barcodes_no_inserts)}\n"
        )

if __name__ == '__main__':
    import sys
    empty_insert_plots(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
