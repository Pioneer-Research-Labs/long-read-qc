process extract_barcodes {
    publishDir "$params.localdir/$meta.id",  mode: 'copy'
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(reads), path(construct), path(flanking)

    output:
    tuple val(meta), path("barcodes.fasta")
    path "cutadapt_barcode_report.json"

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_barcode) \
        --discard-untrimmed \
        --revcomp \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o barcodes_raw.fasta \
        -j $task.cpus \
        --json cutadapt_barcode_report.json \
        $reads
    seqkit seq -j $task.cpus --min-len $params.min_bc_len --max-len $params.max_bc_len \
        barcodes_raw.fasta > barcodes.fasta

    """
}