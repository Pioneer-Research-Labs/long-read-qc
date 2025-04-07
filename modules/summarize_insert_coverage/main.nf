process summarize_insert_coverage{
    publishDir("$params.localdir"),  mode: 'copy'
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Summarizing insert coverage'

    input:
    path insert_coverage_map


    output:
        path ('*.csv', arity: '3')
        path('*.png', arity: '2')


    script:
    """
    summarize_and_plot.py $insert_coverage_map insert_coverage
    """
}