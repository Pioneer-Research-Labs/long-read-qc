process plot_insert_histogram {
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    input:
    tuple val(meta), path(insert_fasta)

    output:
    tuple val(meta), path ('insert_length_histogram.png')

    script:
    """
    seqkit fx2tab -li $insert_fasta > inserts.tsv
    summarize_and_plot.py inserts.tsv insert_histogram
    """
}