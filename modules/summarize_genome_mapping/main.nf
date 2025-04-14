process summarize_genome_mapping {
    publishDir("$params.localdir"),  mode: 'copy'
    publishDir("$params.outdir"),  mode: 'copy'


    cpus params.cores

    input:
    path genome_mapping
    path seq_stats_mapping

    output:
    path "concatenated_genome_mapping.csv"
    path "genome_mapping_summary.csv"

    script:
    """
    summarize_and_plot.py $genome_mapping genome_mapping $seq_stats_mapping
    """
}