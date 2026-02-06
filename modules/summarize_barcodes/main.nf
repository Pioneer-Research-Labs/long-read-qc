process summarize_barcodes {
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Summarizing barcodes'

    input:
    path sample_map

    output:
        path ('*.csv', arity: '1')
        path('*.png', arity: '1')

    script:
    """
    summarize_and_plot.py $sample_map barcode
    """
}