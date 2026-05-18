process empty_insert_histogram {
    publishDir "$params.outdir/$meta.id/summary_and_plots",pattern:"*.csv", mode: 'copy'

    input:
        tuple val(meta), path(barcodes), path(inserts)

    output:
        tuple val("$meta.id"),
        path('counts_of_inserts_barcodes_flanking_site_sequences.csv')

    script:
    """
     seqkit fx2tab -li $barcodes > barcodes.tsv;
     seqkit fx2tab -li $inserts > inserts.tsv;
     empty_insert_plots.py barcodes.tsv inserts.tsv $meta.id
    """
}
