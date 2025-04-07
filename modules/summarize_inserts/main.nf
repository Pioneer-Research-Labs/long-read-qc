process summarize_inserts{
    publishDir("$params.localdir"),  mode: 'copy'
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Summarizing inserts'

    input:
    path insert_map

    output:
        path 'concatenated_inserts.csv'
        path 'insert_length_distribution.csv'
        path 'insert_length_distribution.png'

    script:
    """
    summarize_and_plot.py $insert_map insert
    """
}