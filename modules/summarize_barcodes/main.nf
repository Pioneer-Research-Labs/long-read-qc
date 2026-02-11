process summarize_barcodes {
    publishDir("$params.outdir/summary_and_plots"),  mode: 'copy'
    tag 'Summarizing barcodes'

    input:
    path sample_map

    output:
        path('barcode_copy_number.png')

    script:
    """
    summarize_and_plot.py $sample_map barcode
    """
}