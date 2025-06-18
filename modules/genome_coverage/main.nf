process genome_coverage {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(index), path(stats)

    output:
    tuple val(meta), path('genome_coverage.tsv'), path('genome_cov_stats.tsv')
    script:
    """
    samtools coverage $bam > genome_cov_stats.tsv
    bedtools genomecov -ibam $bam -dz > genome_coverage.tsv
    """

}