process rotate_reads {
    publishDir("$params.localdir/$meta.id"),  mode: 'copy'
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("reads_rotated.fasta")

    script:
    """
    seqkit fq2fa $reads | rotate -s $params.rotate_anchor -m 4 - | \
         seqkit replace -p .+ -r "read_{nr}" > reads_rotated.fasta
    """
}