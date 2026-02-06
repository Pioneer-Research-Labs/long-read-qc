process extract_inserts_with_truncated_flanks {
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(flanking), path(untrimmed_fastq)

    output:
    tuple val(meta), path("inserts_from_truncated_flanks.fasta"),path ("untrimmed_from_truncated_flanks.fastq")

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_insert True) \
        --revcomp \
        --untrimmed-output untrimmed_from_truncated_flanks.fastq \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o inserts_truncated_cutadapt.fasta \
        -j $task.cpus \
        $untrimmed_fastq
    seqkit seq --min-len 1 inserts_truncated_cutadapt.fasta > inserts_from_truncated_flanks.fasta
    """
}