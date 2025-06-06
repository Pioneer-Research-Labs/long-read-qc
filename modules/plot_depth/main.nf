process plot_depth{
    publishDir("$params.outdir/$meta.id") ,  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bam_index),  path(insert_stats)

    output:
    tuple val(meta), path('coverage_plot.png')

    script:
    """
    samtools depth -a $bam > depth_report.tsv
    plot_coverage.py depth_report.tsv coverage_plot.png $meta.id
    """
}