#!/usr/bin/env nextflow


// Run the workflow

workflow {


    input_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            meta = [id:row.id, genome:row.genome, construct:file(row.construct)]
            [meta, file(row.file)]
    
        }



    // --- Prep

    // reorient the reads
    //reads = rotate_reads(input_ch)

    map_vector(input_ch)

    // get the flanking sequences from the .dna file
    flanking = get_flanks(input_ch)
    

    // test cutadapt here
    (barcodes, bc_report, bc_tab) = extract_barcodes(input_ch.join(flanking))
    (inserts, ins_report, in_tab) = extract_inserts(input_ch.join(flanking)) 
    
    //.join(reads)
    input_ch
        .join(inserts)
        .join(barcodes) 
     | seq_stats


    barcode_counts(barcodes)


    // mapping inserts

    // split based on single genome or metagenome mode
    splits = inserts.branch { meta, path ->
        single: meta.genome != 'meta'
            [meta, path]
        multi: meta.genome == 'meta'
            [meta, path]
    }


    mapped = map_inserts(splits.single)
    insert_coverage(mapped)
    
    // metagenomic samples
    sketch(splits.multi) | classify | taxonomy


    // report
    report = channel.fromPath("${projectDir}/assets/report_template.ipynb")
    report_utils = channel.fromPath("${projectDir}/assets/report_utils_template.py")

    prepare_report(report, report_utils)

    channel.fromPath(params.samplesheet) | samples


    
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

process map_vector {
    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('mapped_vector.bam'), path("mapped_vector.bam.bai"), path("mapped_vector_stats.tsv")

    script:
    """
    convert_dna.py $meta.construct | \
        minimap2 -ax map-ont --secondary=no - $reads | samtools view -b - | samtools sort - -o 'mapped_vector.bam'
    samtools index mapped_vector.bam
    samtools flagstat -O tsv mapped_vector.bam > mapped_vector_stats.tsv
    """

}

process seq_stats {
       
    publishDir("$params.outdir/$meta.id")

    tag "$meta.id"

    input:
    tuple val(meta), path(raw), path(ins_seqs), path(bc_seqs)
    
    //path(rotated), 
    
    output:
    path 'seq_stats.tsv'

    script:
    """
    seqkit stats -T $raw $ins_seqs $bc_seqs > seq_stats.tsv
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
        -o inserts_cutadapt.fasta \
        --json cutadapt_inserts_report.json \
        $reads
    seqkit seq --min-len 1 inserts_cutadapt.fasta > inserts.fasta
    seqkit fx2tab -li inserts.fasta > inserts.tsv
    """
}


// Insert mapping

process map_inserts {

    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    input:
    tuple val(meta), path(ins_seqs)

    output:
    tuple val(meta), path('mapped_inserts.bam'), path("mapped_inserts.bam.bai"), path('mapped_insert_stats.tsv')

    script:
    """
    export ref_fa="\$HOME/shared/genomes/${meta.genome}/${meta.genome}_contigs.fna"
    
    minimap2 -ax map-ont \$ref_fa $ins_seqs | samtools view -b - | samtools sort - -o mapped_inserts.bam
    samtools index mapped_inserts.bam
    samtools flagstats -O tsv mapped_inserts.bam > mapped_insert_stats.tsv
    """

}


process insert_coverage {

    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(index), path(stats)

    output:
    tuple val(meta), path('gene_coverage.bed'), path('insert_coverage.bed'), 
        path('genome_coverage.tsv'), path('genome_cov_stats.tsv'), path("insert_coverage_full.bed"),
        path('insert_intersect.out')

    script:
    """
    export gff="\$HOME/shared/genomes/${meta.genome}/${meta.genome}_genes.gff"
    export bed="\$HOME/shared/genomes/${meta.genome}/${meta.genome}_genes.bed"

    bedtools coverage -a \$gff -b $bam > gene_coverage.bed
    bedtools coverage -b \$gff -a <(bedtools bamtobed -i $bam) > insert_coverage.bed
    bedtools coverage -b \$gff -a <(bedtools bamtobed -i $bam) -F 1 > insert_coverage_full.bed
    bedtools intersect -a <(bedtools bamtobed -i $bam) -b \$bed  -wao > insert_intersect.out 
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




// --- Metagenomics

process sketch {
    
    publishDir("$params.outdir/$meta.id")
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

process classify {
    
    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(insert_sig)
  
    output:
    tuple val(meta), path('insert_matches.csv')
  
    script:
    """
    sourmash scripts fastgather -o insert_matches.csv -c $task.cpus -t $params.meta_ovlp -k 21 \
        inserts.sig.gz $params.sourmash_db
    """
}

process taxonomy {
    publishDir("$params.outdir/$meta.id")
    tag "$meta.id"

    input:
    tuple val(meta), path(insert_matches)

    output:
    tuple val(meta), path('insert_taxonomy.csv')

    script:
    """
     sourmash tax metagenome -g $insert_matches -t $params.taxonomy -F csv_summary > insert_taxonomy.csv
    """

}