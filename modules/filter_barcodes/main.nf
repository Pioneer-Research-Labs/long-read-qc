process filter_barcodes {
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path("barcodes_filtered.fasta")

    script:
    """
    seqkit seq -j $task.cpus --min-len $params.min_bc_len --max-len $params.max_bc_len \
        $barcodes > barcodes_filtered.fasta
    """
}