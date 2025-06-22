process process_8bp_genome_tags {
    publishDir("$params.outdir/$meta.id"),  mode: 'move'
    tag "$meta.id"

    input:
    tuple val(meta), path(cutadapt_genome_tags_tsv), path(tesseract_genome_tags), path(fastq), path(construct)

    output:
    path('*.fq.gz')
    path("$meta.id"+"_sample_sheet.csv")

    when:
    script:
    """
    echo $fastq
    process_genome_tags.py $cutadapt_genome_tags_tsv $tesseract_genome_tags $meta.id $fastq $construct $params.outdir
    """

}