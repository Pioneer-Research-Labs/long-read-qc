process get_truncated_inserts_as_tsv {
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    input:
    tuple val(meta), path(truncated_insert_fasta), path (untrimmed_from_truncated_flanks_fq)

    output:
    tuple val(meta), path("truncated_inserts.tsv")

    script:
    """
    seqkit fx2tab -li $truncated_insert_fasta > truncated_inserts.tsv
    """
}