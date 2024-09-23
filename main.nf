#!/usr/bin/env nextflow


// Run the workflow

workflow {


    input_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            meta = [id:row.id, construct:file(row.construct)]
            [meta, file(row.file)]
    
        }



    // --- Prep

    // reorient the reads
    read_ch = rotatereads(input_ch)

    // get the flanking sequences from the .dna file
    flanking_ch = getflanks(input_ch)
    


    // --- find and extract barcodes and inserts
    
    // map flanking sequences
    map_ch = mapflanking(read_ch.join(flanking_ch))

    // convert maps to bed the filter
    coords = extract_coords(map_ch)


    (filtered_ch, stats) = filter_bed(coords)

    
    // use bed to extract sequences
    seqs = extract_seqs(filtered_ch.join(read_ch))

    // convert to tabular
    tab = seqtotab(seqs)
    
    // some stats
    seqstats(seqs)

    // count unique barcodes
    barcodecounts(tab)

    // mapping inserts
    map_ch = mapinserts(seqs, params.ref_fa)
    insertcoverage(map_ch, params.ref_gff)
    
    // report
    channel.fromPath("${projectDir}/assets/report_template.ipynb") \
        | preparereport

}


// Processes

process getflanks {
    publishDir("$params.outdir/$meta.id")
    tag ("Extracting flanks for $meta.id")

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
    tag "Rotating reads for $meta.id"

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
    tuple val(meta), path(reads), path(flanking)

    output:
    tuple val(meta), path('flanking_maps.paf')

    script:
    """
    minimap2 -cxsr -k $params.kmer_size -m $params.chaining_score -B $params.mismatch -N 0 \
        $flanking $reads > flanking_maps.paf 
    """

}

process extract_coords {

    publishDir("$params.outdir/$meta.id")

    tag "Extracting coordinates for $meta.id"

    input:
    tuple val(meta), path(flanking_maps)

    output:
    tuple val(meta), path('barcode_coords.bed'), path('insert_coords.bed')

    script:
    """
    align_coords.py $flanking_maps
    """

}

process filter_bed {
    publishDir("$params.outdir/$meta.id")

    tag "Filtering bed for $meta.id"

    input:
    tuple val(meta), path(barcode_bed), path(insert_bed)

    output:
    tuple val(meta), path('barcode_filtered.bed'), path('insert_filtered.bed')
    tuple val(meta), path('barcode_stats.json'), path('insert_stats.json')

    script:
    """
    filter_bed.py $barcode_bed 'barcode_filtered.bed' 'barcode_stats.json'
    filter_bed.py $insert_bed 'insert_filtered.bed' 'insert_stats.json'
    """
}


process extract_seqs {
    
    publishDir("$params.outdir/$meta.id")

    tag "Extracting sequences for $meta.id"

    input:
    tuple val(meta), path(barcode_bed), path(insert_bed), path(reads)

    output:
    tuple val(meta), path('barcode_seqs.fasta'), path('insert_seqs.fasta')

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
    tuple val(meta), path(bc_seqs), path(ins_seqs)

    output:
    path 'seq_stats.out'

    script:
    """
    seqkit stats -T $bc_seqs $ins_seqs > seq_stats.out
    """
}

process seqtotab {
    
    publishDir("$params.outdir/$meta.id")
    tag 'Sequences to table'

    input:
    tuple val(meta), path(bc_seqs), path(ins_seqs)

    output:
    tuple val(meta), path('barcode_seqs.tsv'), path('insert_seqs.tsv')

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
    tuple val(meta), path(bc_seqs_tab), path(insert_seqs_tab)

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
    tag "Mapping inserts to genome for $meta.id"

    input:
    tuple val(meta), path(bc_seqs), path(ins_seqs)
    path ref

    output:
    tuple val(meta), path('mapped_inserts.bam'), path("mapped_inserts.bam.bai")

    script:
    """
    minimap2 -ax map-ont $ref $ins_seqs | samtools view -b - | samtools sort - -o mapped_inserts.bam
    samtools index mapped_inserts.bam
    """

}


process insertcoverage {

    publishDir("$params.outdir/$meta.id")
    tag "Calculating coverage for $meta.id"

    input:
    tuple val(meta), path(bam), path(index)
    path gff

    output:
    tuple val(meta), path('gene_coverage.bed'), path('insert_coverage.bed'), path('genome_coverage.tsv'), path('genome_cov_stats.tsv')

    script:
    """
    bedtools coverage -a $gff -b $bam > gene_coverage.bed
    bedtools coverage -b $gff -a <(bedtools bamtobed -i $bam) > insert_coverage.bed
    bedtools genomecov -ibam $bam -dz > genome_coverage.tsv
    samtools coverage $bam > genome_cov_stats.tsv
    """

}

process preparereport {

    publishDir("$params.outdir")
    tag 'Preparing report'

    input:
    path report

    output:
    path 'report.ipynb'

    script:
    """
    cp $report 'report.ipynb'
    """
}