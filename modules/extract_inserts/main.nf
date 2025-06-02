process extract_inserts {
    publishDir "$params.localdir/$meta.id",  mode: 'copy'
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(reads), path(construct), path(flanking)

    output:
    tuple val(meta), path("inserts.fasta")
    path ("cutadapt_inserts_report.json")
    path ("cutadapt_info.tsv")
    tuple val(meta), path ("untrimmed.fastq")


    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_insert) \
        --revcomp \
        --untrimmed-output untrimmed.fastq \
        --info-file cutadapt_info.tsv \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o inserts_cutadapt.fasta \
        -j $task.cpus \
        --json cutadapt_inserts_report.json \
        $reads
    seqkit seq --min-len 1 inserts_cutadapt.fasta > inserts.fasta
    """
}