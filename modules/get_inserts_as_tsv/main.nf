process get_inserts_as_tsv {
    tag("$meta.id")

    input:
    tuple val(meta), path(insert_fasta)

    output:
    tuple val(meta), path("inserts.tsv")

    script:
    """
    seqkit fx2tab -li $insert_fasta > inserts.tsv
    """
}