process map_catalog {
    publishDir("$params.outdir/$meta.id/primary_data"), mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(ins_seqs), path(catalog_fasta)

    output:
    tuple val(meta), path('mapped_catalog.bam'), path('mapped_catalog.bam.bai'), path('mapped_catalog_stats.tsv')

    script:
    """
    minimap2 -ax $params.tech -t $task.cpus $catalog_fasta $ins_seqs \
        | samtools view -@ $task.cpus -b - \
        | samtools sort - -@ $task.cpus -o mapped_catalog.bam
    samtools index -@ $task.cpus mapped_catalog.bam
    samtools flagstats -@ $task.cpus -O tsv mapped_catalog.bam > mapped_catalog_stats.tsv
    """
}
