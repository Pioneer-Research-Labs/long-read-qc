process extract_inserts_with_truncated_flanks {
    publishDir "$params.localdir/$meta.id",  mode: 'copy'
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(flanking), path(untrimmed_fastq)

    output:
    tuple val(meta), path("inserts_from_truncated_flanks.fasta"),path ("cutadapt_inserts_report_from_truncated_flanks.json"),
    path ("untrimmed_from_truncated_flanks.fastq"), path ("cutadapt_info_from_truncated_flanks.tsv")

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_insert True) \
        --revcomp \
        --untrimmed-output untrimmed_from_truncated_flanks.fastq \
        --info-file cutadapt_info_from_truncated_flanks.tsv \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o inserts_truncated_cutadapt.fasta \
        -j $task.cpus \
        --json cutadapt_inserts_report_from_truncated_flanks.json \
        $untrimmed_fastq
    seqkit seq --min-len 1 inserts_truncated_cutadapt.fasta > inserts_from_truncated_flanks.fasta
    """
}