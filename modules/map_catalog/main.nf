process map_catalog {
    publishDir("$params.outdir/$meta.id/primary_data"), mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(ins_seqs), path(catalog_fasta)

    output:
    tuple val(meta), path('mapped_inserts.bam'), path('mapped_inserts.bam.bai'), path('mapped_inserts_stats.tsv')

    script:
    """
    minimap2 -ax $params.tech -t $task.cpus $catalog_fasta $ins_seqs \
        | samtools view -@ $task.cpus -b - \
        | samtools sort - -@ $task.cpus -o mapped_inserts.bam
    samtools index -@ $task.cpus mapped_inserts.bam
    samtools flagstats -@ $task.cpus -O tsv mapped_inserts.bam > mapped_inserts_stats.tsv
    """
}
