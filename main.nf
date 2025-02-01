#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
Usage: nextflow run Pioneer-Research-Labs/long-read-qc -latest

Options:
--samplesheet <file>      Path to the sample sheet (default: samplesheet.csv)
--outdir <dir>            Output directory (default: results)
--error_rate <float>      Error rate for barcode searching (default: 0.1)
--min_overlap <int>       Minimum overlap for barcode searching (default: 3)
--min_bc_len <int>        Minimum barcode length (default: 20)
--max_bc_len <int>        Maximum barcode length (default: 60)
--meta_ovlp <int>         Overlap bp for sourmash (default: 1000)
--sourmash_db <file>      Path to the sourmash database (default: /srv/shared/databases/sourmash/gtdb-rs220-k21.zip)
--taxonomy <file>         Path to the taxonomy database (default: /srv/shared/databases/sourmash/gtdb-rs220.lineages.sqldb)
--cores <int>             Number of cores to use (default: 4)

"""
}


// Run the workflow

workflow {

    log.info """
▗▄▄▖▗▄▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▄▖▗▄▄▖     ▗▄▄▖▗▄▄▄▖▗▄▄▖ ▗▄▄▄▖▗▖   ▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖ ▗▄▄▖
▐▌ ▐▌ █  ▐▌ ▐▌▐▛▚▖▐▌▐▌   ▐▌   ▐▌ ▐▌    ▐▌ ▐▌ █  ▐▌ ▐▌▐▌   ▐▌     █  ▐▛▚▖▐▌▐▌   ▐▌   
▐▛▀▘  █  ▐▌ ▐▌▐▌ ▝▜▌▐▛▀▀▘▐▛▀▀▘▐▛▀▚▖    ▐▛▀▘  █  ▐▛▀▘ ▐▛▀▀▘▐▌     █  ▐▌ ▝▜▌▐▛▀▀▘ ▝▀▚▖
▐▌  ▗▄█▄▖▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▙▄▄▖▐▌ ▐▌    ▐▌  ▗▄█▄▖▐▌   ▐▙▄▄▖▐▙▄▄▖▗▄█▄▖▐▌  ▐▌▐▙▄▄▖▗▄▄▞▘

Long Read Processing and QC Pipeline          
"""    

 // Show help message

    if (params.help) {
        helpMessage()
        exit 0
    }

    channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            meta = [id:row.id, genome:row.genome]
            [meta, file(row.file)]
    
        } 
        | set {input_ch}

    channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            meta = [id:row.id, genome:row.genome]
            reads = file(row.file)
            construct = file(row.construct)
            [meta, reads, construct]
        }
        | set {constructs}

    // reorient the reads
    //reads = rotate_reads(input_ch)
    // Generate quality report using fastplong
    quality_report(constructs)
    map_vector(constructs)


    // get the flanking sequences from the .dna file
    flanking = get_flanks(constructs)

    // extract barcodes and inserts
    (barcodes, bc_report, bc_tab) = extract_barcodes(input_ch.join(flanking))
    (inserts, ins_report, in_tab) = extract_inserts(input_ch.join(flanking))

    // combine for read stats
    input_ch
        .join(inserts)
        .join(barcodes)
     | seq_stats


    barcode_counts(barcodes)


    // mapping inserts

    //split based on single genome or metagenome mode
    splits = inserts.branch { meta2, path ->
        single: meta2.genome != 'meta'
            [meta2, path]
        multi: meta2.genome == 'meta'
            [meta2, path]
    }

    mapped = map_inserts(splits.single)
    insert_outputs = insert_coverage(mapped)
    generate_plots(insert_outputs)
    plot_depth(insert_outputs)

    // metagenomic samples
    sketched = sketch(splits.multi)
    classified = classify(sketched, params.sourmash_db)
    taxonomy(classified, params.taxonomy)


    // report
    report = channel.fromPath("${projectDir}/assets/report_template.ipynb")
    report_utils = channel.fromPath("${projectDir}/bin/report_utils_template.py")

    prepare_report(report, report_utils)

    channel.fromPath(params.samplesheet) | samples




}


// Processes

process get_flanks {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag ("$meta.id")

    input:
    tuple val(meta), path(reads), path(construct)

    output:
    tuple val(meta), path("flanking.gb")

    script:
    """
    echo $projectDir
    get_flanking.py $construct
    """
}

process rotate_reads {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
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
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(construct)

    output:
    tuple val(meta), path('mapped_vector.bam'), path("mapped_vector.bam.bai"), path("mapped_vector_stats.tsv")

    script:
    """
    echo $construct
    convert_dna.py $construct | \
        minimap2 -ax map-ont --secondary=no - $reads | samtools view -b - | samtools sort - -o 'mapped_vector.bam'
    samtools index mapped_vector.bam
    samtools flagstat -O tsv mapped_vector.bam > mapped_vector_stats.tsv
    """

}

process seq_stats {
       
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'

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
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
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
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
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

    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
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
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
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

    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    containerOptions '--volume $HOME/shared/genomes:/genomes'

    input:
    tuple val(meta), path(ins_seqs)

    output:
    tuple val(meta), path('mapped_inserts.bam'), path("mapped_inserts.bam.bai"), path('mapped_insert_stats.tsv')

    script:
    """
    export ref_fa="${params.genomes}/${meta.genome}/${meta.genome}_contigs.fna"

    minimap2 -ax map-ont -t $task.cpus \$ref_fa $ins_seqs | samtools view -b - | samtools sort - -o mapped_inserts.bam
    samtools index mapped_inserts.bam
    samtools flagstats -O tsv mapped_inserts.bam > mapped_insert_stats.tsv
    """

}


process insert_coverage {

    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    containerOptions '--volume $HOME/shared/genomes:/genomes'

    input:
    tuple val(meta), path(bam), path(index), path(stats)

    output:
    tuple val(meta), path('gene_coverage.bed'), path('insert_coverage.bed'), 
        path('genome_coverage.tsv'), path('genome_cov_stats.tsv'), path("insert_coverage_full.bed"),
        path('insert_intersect.out'), path('depth_report.tsv')

    script:
    """
    export gff="${params.genomes}/${meta.genome}/${meta.genome}_genes.gff"
    export bed="${params.genomes}/${meta.genome}/${meta.genome}_genes.bed"

    bedtools coverage -a \$gff -b $bam > gene_coverage.bed
    bedtools coverage -b \$gff -a <(bedtools bamtobed -i $bam) > insert_coverage.bed
    bedtools coverage -b \$gff -a <(bedtools bamtobed -i $bam) -F 1 > insert_coverage_full.bed
    bedtools intersect -a <(bedtools bamtobed -i $bam) -b \$bed  -wao > insert_intersect.out 
    bedtools genomecov -ibam $bam -dz > genome_coverage.tsv
    samtools coverage $bam > genome_cov_stats.tsv
    samtools depth -a $bam > depth_report.tsv
    """

}

process prepare_report {

    publishDir("$params.outdir"),  mode: 'copy'
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

    publishDir("$params.outdir"),  mode: 'copy'
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

process classify {
    
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    cpus params.cores

    input:
    tuple val(meta), path(insert_sig)
    path sourmash_db
  
    output:
    tuple val(meta), path('insert_matches.csv')
  
    script:
    """
    sourmash scripts fastgather -o insert_matches.csv -c $task.cpus -t $params.meta_ovlp -k 21 \
        inserts.sig.gz $sourmash_db
    """
}

process taxonomy {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    memory '32 GB'
    input:
    tuple val(meta), path(insert_matches)
    path taxonomy

    output:
    tuple val(meta), path('insert_taxonomy.csv')

    script:
    """
     sourmash tax metagenome -g $insert_matches -t $taxonomy -F csv_summary > insert_taxonomy.csv
    """

}

process quality_report {

    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(construct)

    output:
    tuple val(meta), path('fastplong.html'), path('fastplong.fq')

    script:
    """
    fastplong -i $reads -o fastplong.fq  -A -Q -L
    """
}

process generate_plots {
    publishDir("$params.outdir"),  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path('gene_coverage.bed'), path('insert_coverage.bed'),
        path('genome_coverage.tsv'), path('genome_cov_stats.tsv'), path("insert_coverage_full.bed"),
        path('insert_intersect.out'), path('depth_report.tsv')

    output:
    path('raw_seq_stats.csv')
    path("seq_summary.csv")
    path("barcode_length_distribution.png")
    path("barcode_proportions.png")
    path("barcode_copy_number.png")
    path("insert_length_distribution.png")
    path("full_genes_per_fragment.png")
    path("partial_genes_per_fragment.png")
    path("barcode_lengths.csv")
    path("barcode_proportions.csv")
    path("barcode_copy_number.csv")
    path("insert_length_distribution.csv")
    path("full_genes_per_fragment.csv")
    path("partial_genes_per_frament.csv")
    script:
    """
    echo $PWD
    ls -lt $PWD/results
    visualize_results.py $PWD/$params.outdir/samples.csv  $PWD/$params.outdir
    """
}

process plot_depth{
    publishDir("$params.outdir/$meta.id") ,  mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path('gene_coverage.bed'), path('insert_coverage.bed'),
        path('genome_coverage.tsv'), path('genome_cov_stats.tsv'), path("insert_coverage_full.bed"),
        path('insert_intersect.out'), path('depth_report.tsv')

    output:
    tuple val(meta), path('coverage_plot.png')

    script:
    """
    plot_coverage.py depth_report.tsv coverage_plot.png $meta.id
    """
}
