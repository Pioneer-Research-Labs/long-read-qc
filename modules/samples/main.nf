process samples {
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Moving sample sheet'

    input:
    path samplesheet

    output:
    path 'samples.csv'

    script:
    """
    cp $samplesheet 'samples.csv'
    """
}