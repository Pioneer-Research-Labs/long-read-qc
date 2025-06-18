process map_genome {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(reads), path(construct), path(fna)

    output:
    tuple val(meta), path('mapped_genome.bam'), path("mapped_genome.bam.bai"), path("mapped_genome_stats.tsv")

    script:
    """
    minimap2 -ax $params.tech -t $task.cpus $fna $reads | samtools view -@ $task.cpus -b - | samtools sort - -@ $task.cpus -o 'mapped_genome.bam'
    samtools index -@ $task.cpus mapped_genome.bam
    samtools flagstat -@ $task.cpus -O tsv mapped_genome.bam > mapped_genome_stats.tsv
    """

}