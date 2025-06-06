process quality_report {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(construct)

    output:
    tuple val(meta), path('fastplong.html'), path('fastplong.fq')

    script:
    """
    fastplong -i $reads -o fastplong.fq  -A -Q -L
    """
}