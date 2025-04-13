process summarize_genome_coverage {
    publishDir("$params.localdir"),  mode: 'copy'
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Summarizing genome coverage'

    input:
    path genome_cov_mapping

    output:
        path('genome_coverage.png')


    script:
    """
    summarize_and_plot.py $genome_cov_mapping genome_coverage
    """
}