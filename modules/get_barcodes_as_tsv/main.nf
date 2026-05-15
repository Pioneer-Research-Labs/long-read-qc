process get_barcodes_as_tsv{
    tag("$meta.id")

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path ("barcodes.tsv")

    script:
    """
    seqkit fx2tab -li $barcodes > barcodes.tsv
    """

}