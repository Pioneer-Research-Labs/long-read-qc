process barcode_counts {
    publishDir("$params.outdir/$meta.id/primary_data"),  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val("$meta.id"), path('barcode_counts.tsv')

    script:
    """
     seqkit fx2tab -j $task.cpus -i $barcodes | cut -f2 | sort | uniq -c | \
        awk '{print \$2"\t"\$1}' > barcode_counts.tsv
    """
}