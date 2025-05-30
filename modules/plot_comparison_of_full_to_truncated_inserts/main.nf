process plot_comparison_of_full_to_truncated_inserts {
    publishDir "$params.localdir/$meta.id", mode: 'copy'
    publishDir "$params.outdir/$meta.id", mode: 'copy'
    tag("$meta.id")

    input:
    tuple val(meta), path(full_inserts_fasta), path(truncated_inserts_fasta),
    path(untrimmed_from_truncated_flanks_fastq), path(reads)

    output:
    path("truncated_vs_intact_flanks_comparison_plot.png")
    path("all_seq_ids_lengths.txt")
    path("full_insert_read_lengths.txt")
    path("truncated_ids_read_lengths.txt")
    path("untrimmed_ids_lengths.txt")

    script:
    """
    # Get the seq id and length for all sequences
    seqkit fx2tab -nli $reads  > all_seq_ids_lengths.txt
    # Get the seq id and length for full inserts
    seqkit fx2tab -nli $full_inserts_fasta > full_insert_read_lengths.txt
    # Get the id and length for the truncated inserts
    seqkit fx2tab -nli $truncated_inserts_fasta > truncated_ids_read_lengths.txt
    # Get the seq id and length of all the reads that weren't trimmed with full or truncated flanks
    seqkit fx2tab -nli $untrimmed_from_truncated_flanks_fastq > untrimmed_ids_lengths.txt
    plot_comparison_full_and_truncated_inserts.py full_insert_read_lengths.txt truncated_ids_read_lengths.txt all_seq_ids_lengths.txt untrimmed_ids_lengths.txt ${meta.id}
    """
}