#!/usr/bin/env nextflow


// Run the workflow

workflow {


    input_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            meta = [id:row.id, construct:file(row.construct)]
            [meta, file(row.file)]
    
        }

    input_ch.view()

    ref_ch = channel.fromPath(params.ref_fa)
    gff_ch = channel.fromPath(params.ref_gff)

    // --- Prep

    // reorient the reads
    read_ch = rotatereads(input_ch)

    // get the flanking sequences from the .dna file
    flanking_ch = getflanks(input_ch)
    
    
    // --- find and extract barcodes and inserts

    // map flanking sequences
    map_ch = mapflanking(read_ch, flanking_ch)

    // convert maps to bed
    (bar_coords, ins_coords) = cleanmaps(map_ch)

    // use bed to extract sequences
    (bar_seqs, ins_seqs) = extract_seqs(bar_coords, ins_coords, read_ch)

    // convert to tabular
    (bar_tab, ins_tab) = seqtotab(bar_seqs, ins_seqs)
    
    // some stats
    seqstats(bar_seqs, ins_seqs)

    // count unique barcodes
    barcodecounts(bar_tab)

    // mapping inserts
    map_ch = mapinserts(ins_seqs, ref_ch)
    insertcoverage(map_ch, gff_ch)
    
}


// Processes

process getflanks {
    publishDir("$params.outdir/$meta.id")
    tag ("Extracting flanks")

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("flanking.fasta")

    script:
    """
    get_flanking.py $meta.construct
    """
}

process rotatereads {
    publishDir("$params.outdir/$meta.id")
    tag 'Rotating reads'

    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("reads_rotated.fasta")

    script:
    """
    seqkit fq2fa $reads | rotate -s $params.rotate_anchor -m 4 - > reads_rotated.fasta
    """
}

process mapflanking {
    publishDir("$params.outdir/$meta.id")

    tag "Mapping flanking regions"
    
    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(flanking)

    output:
    tuple val(meta), path('flanking_maps.paf')

    script:
    """
    minimap2 -cxsr -k $params.kmer_size -m $params.chaining_score -B $params.mismatch \
     $flanking $reads > flanking_maps.paf 
    """

}

process cleanmaps {

    publishDir("$params.outdir/$meta.id")

    tag 'Cleaning alignment output'

    input:
    tuple val(meta), path(flanking_maps)

    output:
    tuple val(meta), path('barcode_coords.bed')
    tuple val(meta), path('insert_coords.bed')

    script:
    """
    align_coords.py $flanking_maps
    """

}

process extract_seqs {
    
    publishDir("$params.outdir/$meta.id")

    tag 'Extracting sequences'

    input:
    tuple val(meta), path(barcode_bed)
    tuple val(meta), path(insert_bed)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('barcode_seqs.fasta')
    tuple val(meta), path('insert_seqs.fasta')

    script:
    """
    seqkit subseq --bed $barcode_bed $reads > barcode_seqs.fasta
    seqkit subseq --bed $insert_bed $reads > insert_seqs.fasta
    """

}

process seqstats {
       
    publishDir("$params.outdir/$meta.id")

    tag 'Sequence stats'

    input:
    tuple val(meta), path(bc_seqs)
    tuple val(meta), path(ins_seqs)

    output:
    path 'seq_stats.out'

    script:
    """
    seqkit stats $bc_seqs $ins_seqs > seq_stats.out
    """
}

process seqtotab {
    
    publishDir("$params.outdir/$meta.id")
    tag 'Sequences to table'

    input:
    tuple val(meta), path(bc_seqs)
    tuple val(meta), path(ins_seqs)

    output:
    tuple val(meta), path('barcode_seqs.tsv')
    tuple val(meta), path('insert_seqs.tsv')

    script:
    """
    seqkit fx2tab -l $bc_seqs > 'barcode_seqs.tsv'
    seqkit fx2tab -l $ins_seqs > 'insert_seqs.tsv'
    """
    
}

process barcodecounts {

    publishDir("$params.outdir/$meta.id")
    tag 'Counting unique barcodes'

    input:
    tuple val(meta), path(bc_seqs_tab)

    output:
    path 'barcode_counts.tsv'

    script:
    """
    cut -f 2 $bc_seqs_tab | sort | uniq -c | awk '{print \$2"\t"\$1}' > barcode_counts.tsv
    """
}

// Insert mapping

process mapinserts {

    publishDir("$params.outdir/$meta.id")
    tag 'Mapping inserts to genome'

    input:
    tuple val(meta), path(inserts)
    path ref

    output:
    tuple val(meta), path('mapped_inserts.bam'), path("mapped_inserts.bam.bai")

    script:
    """
    minimap2 -ax map-ont $ref $inserts | samtools view -b - | samtools sort - -o mapped_inserts.bam
    samtools index mapped_inserts.bam
    """

}


process insertcoverage {

    publishDir("$params.outdir/$meta.id")
    tag 'Calculating coverage'

    input:
    tuple val(meta), path(bam), path(index)
    path gff

    output:
    path 'gene_coverage.bed'
    path 'insert_coverage.bed'
    path 'genome_coverage.tsv'

    script:
    """
    bedtools coverage -a $gff -b $bam > gene_coverage.bed
    bedtools coverage -b $gff -a <(bedtools bamtobed -i $bam) > insert_coverage.bed
    bedtools genomecov -ibam $bam -dz > genome_coverage.tsv
    """

}