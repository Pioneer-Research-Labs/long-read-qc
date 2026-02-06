process extract_genome_tags {
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta),  path(construct), path(reads), path(flanking)

    output:
    tuple val(meta), path("genome_tags.fasta"), path ("cutadapt_genome_tags_report.json"),
    path ("cutadapt_genome_tags_info.tsv"), path("extracted_genome_tags.tsv")


    script:
    """

    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_genome_tag False $construct ) \
        --revcomp \
        --info-file cutadapt_genome_tags_info.tsv \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o genome_tags_cutadapt.fasta \
        -j $task.cpus \
        --json cutadapt_genome_tags_report.json \
        $reads
    seqkit seq --min-len 8 --max-len 8 genome_tags_cutadapt.fasta > genome_tags.fasta
    seqkit fx2tab -Q genome_tags.fasta > extracted_genome_tags.tsv
    """
}