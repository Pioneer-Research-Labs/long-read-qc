process empty_insert_histogram {
    publishDir("$params.outdir/$meta.id/summary_and_plots", mode: 'copy')

    input:
        tuple val(meta), path(barcodes), path(inserts), path(sites)

    output:
        tuple val("$meta.id"),
        path("histogram_of_sites_with_barcode_no_insert.png"),
        path('counts_of_inserts_barcodes_flanking_site_sequences.csv'),
        path('sites.tsv')


    script:
    """
     seqkit fx2tab -li $barcodes > barcodes.tsv;
     seqkit fx2tab -li $inserts > inserts.tsv;
     seqkit fx2tab -li $sites > sites.tsv;
     empty_insert_plots.py barcodes.tsv inserts.tsv sites.tsv $meta.id
    """
}