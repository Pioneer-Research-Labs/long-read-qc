#!/usr/bin/env nextflow


// sample sheet
params.samplesheet = "TagLib_nanopore_samples.csv"

// raw nanopore reads
//params.reads = "reads.fastq"

// fasta file with insert and barcode flanking regions
params.flanking = "flanking.fasta"

// output dir
params.outdir = "results"

// vector sequence to use for rotating and orienting reads (using rrnB terminator)
params.rotate_anchor = "TTTATCTGTTGTTTGTCGGTGAACGCTCTC"

// mapping paramaters for finding barcodes and inserts
params.kmer_size = 6
params.chaining_score = 5 
params.mismatch = 4

// genome files for mapping inserts
params.ref_fa = "./data/Helongata.fna"
params.ref_gff = "./data/Helongata.gff"


// Run the workflow

workflow {

    raw_read_ch = channel.fromPath(params.reads)
    flanking_ch = channel.fromPath(params.flanking)
    ref_ch = channel.fromPath(params.ref_fa)
    gff_ch = channel.fromPath(params.ref_gff)
 
    read_ch = rotatereads(raw_read_ch)
    
    // find and extract barcodes and inserts
    map_ch = mapflanking(read_ch, flanking_ch)
    (bar_coords, ins_coords) = cleanmaps(map_ch)
    (bar_seqs, ins_seqs) = extract_seqs(bar_coords, ins_coords, read_ch)
    seqstats(bar_seqs, ins_seqs)
    (bar_tab, ins_tab) = seqtotab(bar_seqs, ins_seqs)
    barcodecounts(bar_tab)

    // mapping inserts
    map_ch = mapinserts(read_ch, ref_ch)
    insertcoverage(map_ch, gff_ch)
}


// Processes

process rotatereads {
    publishDir(params.outdir)
    tag 'Rotating reads'

    input:
    path reads
    
    output:
    path 'reads_rotated.fasta'

    script:
    """
    seqkit fq2fa $reads | rotate -s $params.rotate_anchor -m 4 - > reads_rotated.fasta
    """
}

process mapflanking {
    publishDir(params.outdir)

    tag "Mapping flanking regions"
    
    input:
    path reads
    path flanking

    output:
    path 'flanking_maps.paf'

    script:
    """
    minimap2 -cxsr -k $params.kmer_size -m $params.chaining_score -B $params.mismatch \
     $flanking $reads > flanking_maps.paf 
    """

}

process cleanmaps {

    publishDir(params.outdir)

    tag 'Cleaning alignment output'

    input:
    path flanking_maps

    output:
    path 'barcode_coords.bed'
    path 'insert_coords.bed'

    script:
    """
    align_coords.py $flanking_maps
    """

}

process extract_seqs {
    
    publishDir(params.outdir)

    tag 'Extracting sequences'

    input:
    path barcode_bed
    path insert_bed
    path reads

    output:
    path 'barcode_seqs.fasta'
    path 'insert_seqs.fasta'

    script:
    """
    seqkit subseq --bed $barcode_bed $reads > barcode_seqs.fasta
    seqkit subseq --bed $insert_bed $reads > insert_seqs.fasta
    """

}

process seqstats {
       
    publishDir(params.outdir)

    tag 'Sequence stats'

    input:
    path bc_seqs
    path ins_seqs

    output:
    path 'seq_stats.out'

    script:
    """
    seqkit stats $bc_seqs $ins_seqs > seq_stats.out
    """
}

process seqtotab {
    
    publishDir(params.outdir)   
    tag 'Sequences to table'

    input:
    path bc_seqs
    path ins_seqs

    output:
    path 'barcode_seqs.tsv'
    path 'insert_seqs.tsv'

    script:
    """
    seqkit fx2tab -l $bc_seqs > 'barcode_seqs.tsv'
    seqkit fx2tab -l $ins_seqs > 'insert_seqs.tsv'
    """
    
}

process barcodecounts {

    publishDir(params.outdir)
    tag 'Counting unique barcodes'

    input:
    path bc_seqs_tab

    output:
    path 'barcode_counts.tsv'

    script:
    """
    cut -f 2 $bc_seqs_tab | sort | uniq -c | awk '{print \$2"\t"\$1}' > barcode_counts.tsv
    """
}

// Insert mapping

process mapinserts {

    publishDir(params.outdir)
    tag 'Mapping inserts to genome'

    input:
    path inserts
    path ref

    output:
    path 'mapped_inserts.bam'

    script:
    """
    minimap2 -ax map-ont $ref $inserts | samtools view -b - | samtools sort - -o mapped_inserts.bam
    """

}

process insertcoverage {

    publishDir(params.outdir)
    tag 'Calculating coverage'

    input:
    path bam
    path gff

    output:
    path 'genome_coverage.bed'
    path 'genome_coverage.tsv'
    path 'gene_stats.tsv'

    script:
    """
    bamtocov -r $gff -t 'CDS' --report gene_stats.tsv $bam > genome_coverage.bed
    bedtools genomecov -ibam $bam -dz > genome_coverage.tsv
    """

}