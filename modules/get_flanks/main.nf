process get_flanks {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag ("$meta.id")

    input:
    tuple val(meta), path(construct)

    output:
    tuple val(meta), path("flanking.gb")

    script:
    """
    get_flanking.py $construct
    """
}