process sketch {

    publishDir("$params.localdir/$meta.id"),  mode: 'copy'
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(ins_seqs)

    output:
    tuple val(meta), path('inserts.sig.gz')

    script:
    """
    sourmash sketch dna -p k=21,abund $ins_seqs -o inserts.sig.gz --name inserts
    """
}