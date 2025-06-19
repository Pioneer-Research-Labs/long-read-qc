process build_genome_tags {
    publishDir("$params.outdir/$meta.id"), mode: 'copy'
    tag ("$meta.id")

    input:
    tuple val(meta), path(construct), path(tesseract_table)

    output:
    tuple val(meta), path("combined_tags.csv"), path("prefixed_tags.csv")

    script:
    """
    build_genome_tags.py $construct $tesseract_table
    """
}