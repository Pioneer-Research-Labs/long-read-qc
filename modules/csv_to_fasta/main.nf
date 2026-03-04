process csv_to_fasta {
    tag "$meta.id"

    input:
    tuple val(meta), path(library_csv)

    output:
    tuple val(meta), path('library.fasta')

    script:
    """
    csv_to_fasta.py $library_csv library.fasta
    """
}
