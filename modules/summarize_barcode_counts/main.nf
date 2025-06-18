process summarize_barcode_counts{
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Summarizing barcode counts'

    input:
    path sample_map


    output:
    path 'concatenated_barcode_counts.csv'

    script:
    """
    summarize_and_plot.py $sample_map barcode_counts
    """
}