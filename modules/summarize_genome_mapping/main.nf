process summarize_genome_mapping {
    publishDir("$params.outdir/summary_and_plots"),  mode: 'copy'

    cpus params.cores

    input:
    path genome_mapping
    path seq_stats_mapping

    output:
    path "genome_mapping_summary.csv"

    script:
    """
    summarize_and_plot.py $genome_mapping genome_mapping $seq_stats_mapping
    """
}