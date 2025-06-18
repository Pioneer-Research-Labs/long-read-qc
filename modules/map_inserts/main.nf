process map_inserts {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(ins_seqs), path(fna)

    output:
    tuple val(meta), path('mapped_inserts.bam'), path("mapped_inserts.bam.bai"), path('mapped_insert_stats.tsv')

    script:
    """

    minimap2 -ax $params.tech -t $task.cpus $fna $ins_seqs | samtools view -@ $task.cpus -b - | samtools sort - -@ $task.cpus -o mapped_inserts.bam
    samtools index -@ $task.cpus mapped_inserts.bam
    samtools flagstats -@ $task.cpus -O tsv mapped_inserts.bam > mapped_insert_stats.tsv
    """

}