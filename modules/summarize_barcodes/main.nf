process summarize_barcodes {
    publishDir("$params.localdir"),  mode: 'copy'
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Summarizing barcodes'

    input:
    path sample_map


    output:
        path ('*.csv', arity: '4')
        path('*.png', arity: '3')


    script:
    """
    summarize_and_plot.py $sample_map barcode
    """
}