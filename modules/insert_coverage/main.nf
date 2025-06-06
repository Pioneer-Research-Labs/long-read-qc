process insert_coverage {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(index), path(stats), path(gff), path(bed)

    output:
    tuple val(meta), path('gene_coverage.bed'), path('insert_coverage.bed'), path("insert_coverage_full.bed"),
        path('insert_intersect.out'), path('depth_report.tsv')
    script:
    """

    bedtools coverage -a $gff -b $bam > gene_coverage.bed
    bedtools coverage -b $gff -a <(bedtools bamtobed -i $bam) > insert_coverage.bed
    bedtools coverage -b $gff -a <(bedtools bamtobed -i $bam) -F 1 > insert_coverage_full.bed
    bedtools intersect -a <(bedtools bamtobed -i $bam) -b $bed  -wao > insert_intersect.out
    samtools depth -a $bam > depth_report.tsv
    """

}