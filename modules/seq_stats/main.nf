process seq_stats {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'

    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(raw), path(construct), path(ins_seqs), path(bc_seqs)

    //path(rotated),

    output:
    tuple val("$meta.id"), path ('seq_stats.tsv')

    script:
    """
    seqkit stats -j $task.cpus -T $raw $ins_seqs $bc_seqs > seq_stats.tsv
    """
}