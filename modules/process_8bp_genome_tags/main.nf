process process_8bp_genome_tags {
    publishDir "$params.outdir/$meta.id/summary_and_plots",pattern:"*.tsv", mode: 'copy'
    publishDir "$params.outdir/$meta.id/genome_fastqs",pattern:"*.fq.gz", mode: 'copy'
    publishDir "$params.outdir/$meta.id/summary_and_plots", pattern: "*.png", mode: 'copy'
    publishDir "$params.outdir/$meta.id/primary_data", pattern: "*.csv", mode: 'copy'

    tag "$meta.id"

    input:
    tuple val(meta), path(cutadapt_genome_tags_tsv), path(tesseract_genome_tags), path(fastq), path(construct)

    output:
    tuple val(meta), path('*.fq.gz'), path("$meta.id"+"_sample_sheet.csv"), path("combined_seq_summary.tsv"), path("num_sequences_per_genome.png")

    when:
    script:
    """
    process_genome_tags.py $cutadapt_genome_tags_tsv $tesseract_genome_tags $meta.id $fastq $construct $params.outdir
    """

}