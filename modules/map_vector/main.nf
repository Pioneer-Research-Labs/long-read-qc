process map_vector {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(reads), path(construct)

    output:
    tuple val(meta), path('mapped_vector.bam'), path("mapped_vector.bam.bai"), path("mapped_vector_stats.tsv")

    script:
    """
    echo $construct
    convert_dna.py $construct | \
        minimap2 -ax $params.tech -t $task.cpus --secondary=no - $reads | samtools view -@ $task.cpus -b - | samtools sort - -@ $task.cpus -o 'mapped_vector.bam'
    samtools index -@ $task.cpus mapped_vector.bam
    samtools flagstat -@ $task.cpus -O tsv mapped_vector.bam > mapped_vector_stats.tsv
    """

}