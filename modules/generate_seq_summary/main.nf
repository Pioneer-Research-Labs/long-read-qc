process generate_seq_summary{
    publishDir("$params.localdir"),  mode: 'copy'
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Summarizing sequence stats'

    input:
    path seq_stats_map
    path barcode_map
    path vector_map
    path insert_map

    output:
        path 'seq_summary.csv'
        path 'concatenated_seq_stats.csv'
        path 'concatenated_vector_map_stats.csv'

    script:
    """
    summarize_and_plot.py $seq_stats_map seq_stat $barcode_map $vector_map $insert_map
    """
}