process get_sites_as_tsv {
    publishDir "$params.outdir/$meta.id/primary_data",  mode: 'copy'
    tag("$meta.id")

    input:
    tuple val(meta), path(sites_fasta)

    output:
    tuple val(meta), path("sites.tsv")

    script:
    """
    seqkit fx2tab -li $sites_fasta > sites.tsv
    """
}
