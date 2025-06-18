process taxonomy {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    memory '32 GB'
    input:
    tuple val(meta), path(insert_matches)
    path taxonomy

    output:
    tuple val(meta), path('insert_taxonomy.csv')

    script:
    """
     sourmash tax metagenome -g $insert_matches -t $taxonomy -F csv_summary > insert_taxonomy.csv
    """

}