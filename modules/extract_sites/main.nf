process extract_sites {
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores


    input:
    tuple val(meta), path(reads), path(construct), path(flanking)

    output:
    tuple val(meta), path("sites.fasta")
    path ("extract_sites_report.json")
    path ("extract_sites_info.tsv")
    tuple val(meta), path ("untrimmed_sites.fastq")

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_site) \
        --revcomp \
        --untrimmed-output untrimmed_sites.fastq \
        --info-file extract_sites_info.tsv \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o sites_cutadapt.fasta \
        -j $task.cpus \
        --json extract_sites_report.json \
        $reads
    seqkit seq --min-len 0 sites_cutadapt.fasta > sites.fasta
    """
}