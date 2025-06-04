#!/usr/bin/env nextflow


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES  / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { quality_report } from './modules/quality_report'
include { map_vector } from './modules/map_vector'
include { seq_stats } from './modules/seq_stats'
include { extract_barcodes } from './modules/extract_barcodes'
include { extract_inserts } from './modules/extract_inserts'
include { extract_inserts_with_truncated_flanks } from './modules/extract_inserts_with_truncated_flanks'
include { map_inserts } from './modules/map_inserts'
include { genome_coverage } from './modules/genome_coverage'
include { insert_coverage } from './modules/insert_coverage'
include { map_genome } from './modules/map_genome'
include { sketch } from './modules/sketch'
include { classify } from './modules/classify'
include { taxonomy } from './modules/taxonomy'
include { summarize_barcode_counts } from './modules/summarize_barcode_counts'
include { summarize_inserts } from './modules/summarize_inserts'
include { summarize_barcodes } from './modules/summarize_barcodes'
include { summarize_insert_coverage } from './modules/summarize_insert_coverage'
include { summarize_genome_coverage } from './modules/summarize_genome_coverage'
include { summarize_genome_mapping } from './modules/summarize_genome_mapping'
include { generate_seq_summary } from './modules/generate_seq_summary'
include { plot_comparison_of_full_to_truncated_inserts } from './modules/plot_comparison_of_full_to_truncated_inserts'
include { plot_depth } from './modules/plot_depth'
include { get_flanks } from './modules/get_flanks'
include { get_barcodes_as_tsv } from './modules/get_barcodes_as_tsv'
include { get_inserts_as_tsv } from './modules/get_inserts_as_tsv'
include { get_truncated_inserts_as_tsv } from './modules/get_truncated_inserts_as_tsv'
include { barcode_counts } from './modules/barcode_counts'
include { prepare_report } from './modules/prepare_report'
include { samples } from './modules/samples'


#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    output:
        path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}

def helpMessage() {
    log.info """
Usage: nextflow run Pioneer-Research-Labs/long-read-qc -latest

Options:
--samplesheet <file>      Path to the sample sheet (default: samplesheet.csv)
--outdir <dir>            Output directory (default: results)
--tech <str>              Sequencing technology, map-ont/map-pb/map-hifi (default: map-ont)
--map_genome <bool>       Map all reads to genome (default: false)
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

//  // Show help message
//
//     if (params.help) {
//         helpMessage()
//         exit 0
//     }
//
//     channel.fromPath(params.samplesheet)
//         .splitCsv(header:true)
//         .map { row ->
//             meta = [id:row.id, genome:row.genome]
//             [meta, file(row.file), file(params.constructs + row.construct)]
//
//         }
//         | set {input_ch}
//
//     channel.fromPath(params.samplesheet)
//         .splitCsv(header:true)
//         .map { row ->
//             meta = [id:row.id, genome:row.genome]
//             [meta, file(params.constructs + row.construct)]
//
//         }
//         | set {constructs}
//
//     // Channel holding the sample name, genome and original read paths
//     input_ch
//         .map { meta, reads, construct ->
//             // Return the meta data, reads, and the construct file
//             [meta, reads]
//         } | set {read_ch}
//
//     //Generate quality report using fastplong
//     quality_report(input_ch)
//
//     // This type of pattern collects the results of all the sequences being processed, creates a temp file
//     // and then returns the path to the temp file. This is useful for collecting results from multiple processes
//     // and we do so here since we are collecting the results of the mapping process in a summary table and/or plot.
//     // The vector map is a file that maps the sample to the path of the mapped vector results. The temp files
//     // generated by this method are found in the work/tmp directory.
//     vector_map = map_vector(input_ch).collectFile(){
//         meta, bam, bai, stats ->
//         ["mapped_vector_map.tsv", "${meta.id}\t${params.path_prefix}${stats}\n"]
//     }
//
//
//     // get the flanking sequences from the .dna file
//     flanking = get_flanks(constructs)
//
//     joinChannel = input_ch.join(flanking)
//
//     // extract barcodes
//     (barcodes, bc_report, bc_tab) = extract_barcodes(joinChannel)
//     // extract inserts, returning insert_fasta (with metadata), cutadapt report, cutadapt info, and a fastq file of reads that weren't trimmed
//     (inserts, ins_report, in_tab, untrimmed_meta) = extract_inserts(joinChannel)
//
//
//     // Join untrimmed_meta with input_ch, yielding a channel that contains the metadata,
//     // reads, construct, flanking sequence and the fastq file of reads that weren't trimmed.
//     joinChannel.join(untrimmed_meta)
//         .map { meta, reads, construct, flanking, untrimmed_fq ->
//             // Return the meta data, flanking sequence, and the fastq file of reads that weren't trimmed
//             [meta, flanking, untrimmed_fq]
//         } | set {untrimmed_with_reads}
//
//     // extract any inserts with truncated flanking sequences
//     truncated_data = extract_inserts_with_truncated_flanks(untrimmed_with_reads)
//
//     // Generate a .tsv of the truncated insert data
//     get_truncated_inserts_as_tsv(truncated_data)
//
//     // We need the insert fasta file for the full inserts and the truncated inserts to compare them.
//     data_for_plot = truncated_data.join(inserts)
//          .map { meta, truncated_insert_fasta, cutadapt_truncated_report,  untrimmed_from_truncated_fq, cutadapt_truncated_info, full_inserts ->
//              // Return the meta data, full inserts, truncated inserts, and untrimmed fastq file
//              [meta, full_inserts, truncated_insert_fasta, untrimmed_from_truncated_fq]
//          } // We also need the original sequence file to generate a list of sequence ids/lengths.
//          .join(read_ch)
//
//     // Plot the comparisons between full and truncated inserts
//    plot_comparison_of_full_to_truncated_inserts(data_for_plot)
//
//
//     // combine for read stats
//     combined_data = input_ch
//         .join(inserts)
//         .join(barcodes)
//
//
//     // The seq stats map is a file that maps the sample to the path of the seq stats results file for each sample.
//     seq_stats_results = seq_stats(combined_data).collectFile(){
//         id, tsv ->
//         ["seq_stats_map.tsv", "${id}\t${params.path_prefix}${tsv}\n"]
//         }
//
//     // The sample map is a file that maps the sample to the path of the barcode count file for each sample.
//     barcode_count_sample_map = barcode_counts(barcodes).collectFile(){ id, tsv ->
//         ["sample_map.tsv", "${id}\t${params.path_prefix}${tsv}\n"]
//     }
//
//     // The barcode map is a file that maps the sample to the path of the barcode file for each sample.
//     barcode_map =  get_barcodes_as_tsv(barcodes).collectFile(){ meta, fasta ->
//         ["barcode_map.tsv", "${meta.id}\t${params.path_prefix}${fasta}\n"]
//     }
//
//     // The insert map is a file that maps the sample to the path of the insert file for each sample.
//     insert_map =  get_inserts_as_tsv(inserts).collectFile(){ meta, fasta ->
//         ["insert_map.tsv", "${meta.id}\t${params.path_prefix}${fasta}\n"]
//     }
//
//     //split based on single genome or metagenome mode
//     splits = inserts.branch { meta2, path ->
//         single: meta2.genome != 'meta'
//             [meta2, path]
//         multi: meta2.genome == 'meta'
//             [meta2, path]
//     }
//
//     // map inserts and add the dynamically generated path to the contigs.fna file, adding it to the channel
//     mapped = map_inserts(splits.single | map {
// 	meta, seq_path -> [meta, seq_path, "${params.genomes}/${meta.genome}/${meta.genome}_contigs.fna".toString()]
// 	})
//
//     // Map genome coverage
//     genome_cov_map = genome_coverage(mapped).collectFile(){
//         meta, cov, stats ->
//         ["genome_coverage.tsv", "${meta.id}\t${params.path_prefix}${stats}\n"]
//     }
//     // Add the files needed for insert coverage by dynamically creating paths based on genome.
//     mapped_with_references = mapped | map {
// 	meta, bam, bai, stats -> [
// 		meta, bam,  bai, stats , "${params.genomes}/${meta.genome}/${meta.genome}_genes.gff".toString(),
// 		"${params.genomes}/${meta.genome}/${meta.genome}_genes.bed".toString()]
// 	}
//
//     // get insert coverage and collect the results into a file that maps the sample to the location of
//     // the insert coverage output.
//     insert_outputs = insert_coverage(mapped_with_references
// 			).collectFile(){
//             meta, gene_cov, insert_cov  , insert_cov_full, insert_intersect, depth ->
//         ["insert_coverage.tsv", "${meta.id}\t${params.path_prefix}${insert_cov_full}\n"]
//     }
//
//      // Here we generate various summary files and plots for all the sequences processed
//      summarize_barcode_counts(barcode_count_sample_map)
//      summarize_genome_coverage(genome_cov_map)
//      summarize_inserts(insert_map)
//      summarize_barcodes(barcode_map)
//      seq_summary_results = generate_seq_summary(seq_stats_results, barcode_map, vector_map, insert_map)
//      summarize_insert_coverage(insert_outputs)
//
//     // If we want to, map all reads to the donor genome
//     // Add the genome .fna file to the input ch
//     if (params.map_genome) {
//         genome_map = map_genome(input_ch | map {
//             meta, reads, construct -> [meta, reads, construct, "${params.genomes}/${meta.genome}/${meta.genome}_contigs.fna".toString()]
//         }).collectFile(){
//             meta, bam, bai, stats ->
//             ["mapped_genome_map.tsv", "${meta.id}\t${params.path_prefix}${stats}\n"]
//         }
//         summarize_genome_mapping(genome_map, seq_summary_results[1])
//     }
//      // Plot coverage depth
//      plot_depth(mapped)
//
//     // metagenomic samples
//     sketched = sketch(splits.multi)
//     classified = classify(sketched, params.sourmash_db)
//     taxonomy(classified, params.taxonomy)
//
//
//     // report
//     report = channel.fromPath("${projectDir}/assets/report_template.ipynb")
//     report_utils = channel.fromPath("${projectDir}/assets/report_utils_template.py")
//
//     prepare_report(report, report_utils)
//
//     channel.fromPath(params.samplesheet) | samples
    sayHello()
}
