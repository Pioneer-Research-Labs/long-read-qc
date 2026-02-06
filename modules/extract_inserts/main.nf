process extract_inserts {
    publishDir "$params.outdir/$meta.id/primary_data",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(reads), path(construct), path(flanking)

    output:
    tuple val(meta), path("inserts.fasta")
    tuple val(meta), path ("untrimmed.fastq")


    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_insert ) \
        --revcomp \
        --untrimmed-output untrimmed.fastq \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o inserts_cutadapt.fasta \
        -j $task.cpus  \
        $reads
    seqkit seq --min-len 0 inserts_cutadapt.fasta > inserts.fasta
    """
}