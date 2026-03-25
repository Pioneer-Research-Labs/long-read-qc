process catalog_summary {
    publishDir("$params.outdir/$meta.id/summary_and_plots"), mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai), path(catalog_csv)

    output:
    tuple val(meta), path('catalog_coverage_summary.tsv')

    script:
    """
    samtools idxstats $bam > idxstats.tsv
    catalog_summary.py $catalog_csv idxstats.tsv catalog_coverage_summary.tsv
    """
}
