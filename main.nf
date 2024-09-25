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
    reads = rotate_reads(input_ch)

    seq_stats(reads)

    // get the flanking sequences from the .dna file
    flanking = get_flanks(reads)
    

    // test cutadapt here
    (barcodes, bc_report, bc_tab) = extract_barcodes(reads.join(flanking))
    (inserts, ins_report, in_tab) = extract_inserts(reads.join(flanking)) 


    barcode_counts(barcodes)

    // mapping inserts
    mapped = map_inserts(inserts, params.ref_fa)
    insert_coverage(mapped, params.ref_gff, params.ref_bed)
    
    // report
    report = channel.fromPath("${projectDir}/assets/report_template.ipynb")
    report_utils = channel.fromPath("${projectDir}/assets/report_utils_template.py")

    prepare_report(report, report_utils)

    channel.fromPath(params.samplesheet) | samples


    // --- find and extract barcodes and inserts
    
    /*
    // map flanking sequences
    map_ch = mapflanking(reads.join(flanking))

    // convert maps to bed the filter
    coords = extract_coords(map_ch)


    (filtered_ch, stats) = filter_bed(coords)

    
    // use bed to extract sequences
    seqs = extract_seqs(filtered_ch.join(reads))

    // convert to tabular
    tab = seqtotab(seqs)
    
    // some stats
    seqstats(seqs)

    // count unique barcodes
    barcodecounts(tab)
    */

    
}


// Processes

process get_flanks {
    publishDir("$params.outdir/$meta.id")
    tag ("$meta.id")

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("flanking.fasta")

    script:
    """
    get_flanking.py $meta.construct
    """
}

process rotate_reads {
    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("reads_rotated.fasta")

    script:
    """
    seqkit fq2fa $reads | rotate -s $params.rotate_anchor -m 4 - | \
         seqkit replace -p .+ -r "read_{nr}" > reads_rotated.fasta
    """
}

process seq_stats {
       
    publishDir("$params.outdir/$meta.id")

    tag "$meta.id"

    input:
    tuple val(meta), path(seqs)

    output:
    path 'seq_stats.tsv'

    script:
    """
    seqkit stats -T $seqs > seq_stats.tsv
    """
}

process extract_barcodes {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(reads), path(flanking)

    output:
    tuple val(meta), path("barcodes.fasta")
    path "cutadapt_barcode_report.json"
    path "barcodes.tsv"

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_barcode) \
        --discard-untrimmed \
        --revcomp \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o barcodes_raw.fasta \
        --json cutadapt_barcode_report.json \
        $reads
    seqkit seq --min-len $params.min_bc_len --max-len $params.max_bc_len \
        barcodes_raw.fasta > barcodes.fasta
    seqkit fx2tab -li barcodes.fasta > barcodes.tsv
    """
}

process filter_barcodes {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path("barcodes_filtered.fasta")

    script:
    """
    seqkit seq --min-len $params.min_bc_len --max-len $params.max_bc_len \
        $barcodes > barcodes_filtered.fasta
    """
}

process barcode_counts {

    publishDir("$params.outdir/$meta.id")
    tag("$meta.id")

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path('barcode_counts.tsv')

    script:
    """
     seqkit fx2tab -i $barcodes | cut -f2 | sort | uniq -c | \
        awk '{print \$2"\t"\$1}' > barcode_counts.tsv
    """
}

process extract_inserts {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(reads), path(flanking)

    output:
    tuple val(meta), path("inserts.fasta")
    path "cutadapt_inserts_report.json"
    path "inserts.tsv"

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt_insert) \
        --discard-untrimmed \
        --revcomp \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o inserts.fasta \
        --json cutadapt_inserts_report.json \
        $reads
    seqkit fx2tab -li inserts.fasta > inserts.tsv
    """
}


// Insert mapping

process map_inserts {

    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    input:
    tuple val(meta), path(ins_seqs)
    path ref

    output:
    tuple val(meta), path('mapped_inserts.bam'), path("mapped_inserts.bam.bai")

    script:
    """
    minimap2 -ax map-ont $ref $ins_seqs | samtools view -b - | samtools sort - -o mapped_inserts.bam
    samtools index mapped_inserts.bam
    """

}


process insert_coverage {

    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(index)
    path gff
    path bed

    output:
    tuple val(meta), path('gene_coverage.bed'), path('insert_coverage.bed'), 
        path('genome_coverage.tsv'), path('genome_cov_stats.tsv'), path("insert_coverage_full.bed"),
        path('insert_intersect.out')

    script:
    """
    bedtools coverage -a $gff -b $bam > gene_coverage.bed
    bedtools coverage -b $gff -a <(bedtools bamtobed -i $bam) > insert_coverage.bed
    bedtools coverage -b $gff -a <(bedtools bamtobed -i $bam) -F 1 > insert_coverage_full.bed
    bedtools intersect -a <(bedtools bamtobed -i $bam) -b $bed  -wao > insert_intersect.out 
    bedtools genomecov -ibam $bam -dz > genome_coverage.tsv
    samtools coverage $bam > genome_cov_stats.tsv
    """

}

process prepare_report {

    publishDir("$params.outdir")
    tag 'Preparing report'

    input:
    path report
    path report_utils

    output:
    path 'report.ipynb'
    path 'report_utils.py'

    script:
    """
    cp $report 'report.ipynb'
    cp $report_utils 'report_utils.py'
    """
}

process samples {

    publishDir("$params.outdir")
    tag 'Moving sample sheet'

    input:
    path samplesheet

    output:
    path 'samples.csv'

    script:
    """
    cp $samplesheet 'samples.csv'
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

