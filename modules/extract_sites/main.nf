process extract_sites {
    publishDir "$params.outdir/$meta.id/primary_data",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores


    input:
    tuple val(meta), path(reads), path(construct), path(flanking)

    output:
    tuple val(meta), path("sites.fasta")

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_site) \
        --revcomp \
        -e $params.error_rate \
        -O $params.min_overlap \
        --discard-untrimmed \
        -o sites_cutadapt.fasta \
        -j $task.cpus \
        $reads
    seqkit seq --min-len 0 sites_cutadapt.fasta > sites.fasta
    """
}