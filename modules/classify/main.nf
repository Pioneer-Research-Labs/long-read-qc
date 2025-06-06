process classify {

    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(insert_sig)
    path sourmash_db

    output:
    tuple val(meta), path('insert_matches.csv')

    script:
    """
    sourmash scripts fastgather -o insert_matches.csv -c $task.cpus -t $params.meta_ovlp -k 21 \
        inserts.sig.gz $sourmash_db
    """
}