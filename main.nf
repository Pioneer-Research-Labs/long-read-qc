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
include { extract_sites } from './modules/extract_sites'
include { empty_insert_histogram } from './modules/plot_empty_insert_histogram'
include { extract_inserts_with_truncated_flanks } from './modules/extract_inserts_with_truncated_flanks'
include { extract_genome_tags } from './modules/extract_genome_tags'
include { map_inserts } from './modules/map_inserts'
include { map_genome } from './modules/map_genome'
include { summarize_barcode_counts } from './modules/summarize_barcode_counts'
include { summarize_inserts } from './modules/summarize_inserts'
include { summarize_barcodes } from './modules/summarize_barcodes'
include { summarize_genome_mapping } from './modules/summarize_genome_mapping'
include { generate_seq_summary } from './modules/generate_seq_summary'
include { plot_comparison_of_full_to_truncated_inserts } from './modules/plot_comparison_of_full_to_truncated_inserts'
include { plot_depth } from './modules/plot_depth'
include { get_flanks } from './modules/get_flanks'
include { get_barcodes_as_tsv } from './modules/get_barcodes_as_tsv'
include { plot_insert_histogram } from './modules/plot_insert_histogram'
include { get_inserts_as_tsv } from './modules/get_inserts_as_tsv'
include { get_truncated_inserts_as_tsv } from './modules/get_truncated_inserts_as_tsv'
include { barcode_counts } from './modules/barcode_counts'
include { process_8bp_genome_tags} from './modules/process_8bp_genome_tags'
include { samplesheetToList } from 'plugin/nf-schema'

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
--cores <int>             Number of cores to use (default: 4)
--preprocess_genome_tags  <bool>  Preprocess genome tags (default: false)
"""
}


// Named workflow for pre-processing genome tags
workflow preprocess_genome_tags {

    // Validate the samplesheet
    data = Channel.fromList(samplesheetToList(params.tesseract_samplesheet, "assets/tesseract_samplesheet_validation_schema.json"))
        .map { data ->
            meta = data[0]
            id = [id:meta.id]
            meta = [id, file(params.constructs + meta.construct),file(meta.fastq, checkIfExists:true)]

    }

    genome_constructs = data
        .map { data ->
            // Return the meta data and the construct file
            [data[0], data[1]]

        }

    genome_flanking = get_flanks(genome_constructs)

    joinChannel = data.join(genome_flanking)

    tagChannel = extract_genome_tags(joinChannel)
    process_channel = tagChannel.join(data)

    process_inputs = process_channel
        .map{
            meta, genome_tags, report, genome_tags_info, genome_tags_tsv, construct, fastq ->
            [meta, genome_tags_tsv, file(params.tesseract_oligo_file), fastq, construct]
        }
    sample_sheet_map = process_8bp_genome_tags(process_inputs).collectFile(){
        meta, genome_fastqs, sample_sheet ->
        ["sample_sheet_map.tsv", "${meta.id}\t${params.path_prefix}${sample_sheet}\n"]
    }

    final_sample_sheet = aggregate_sample_sheets(sample_sheet_map)

}

workflow long_read_qc{

        input_ch = Channel.fromList(samplesheetToList(params.samplesheet, "assets/samplesheet_validation_schema.json"))
            .map { meta, construct, sequence ->
                file(params.genomes + meta.genome + "/*_genomic.fna", checkIfExists:true)
                [meta, file(sequence), file(params.constructs + construct, checkIfExists:true)]

            }

        channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row ->
            meta = [id:row.id, genome:row.genome]
            [meta, file(params.constructs + row.construct)]

        }
        | set {constructs}

    // Channel holding the sample name, genome and original read paths
    input_ch
        .map { meta, reads, construct ->
            // Return the meta data, reads, and the construct file
            [meta, reads]
        } | set {read_ch}

    //Generate quality report using fastplong
    quality_report(input_ch)

    // This type of pattern collects the results of all the sequences being processed, creates a temp file
    // and then returns the path to the temp file. This is useful for collecting results from multiple processes
    // and we do so here since we are collecting the results of the mapping process in a summary table and/or plot.
    // The vector map is a file that maps the sample to the path of the mapped vector results. The temp files
    // generated by this method are found in the work/tmp directory.
    vector_map = map_vector(input_ch).collectFile(){
        meta, bam, bai, stats ->
        ["mapped_vector_map.tsv", "${meta.id}\t${params.path_prefix}${stats}\n"]
    }


    // get the flanking sequences from the .dna file
    flanking = get_flanks(constructs)

    joinChannel = input_ch.join(flanking)

    //extract barcodes
    (barcodes, bc_report, bc_tab) = extract_barcodes(joinChannel)

    // extract inserts, returning insert_fasta (with metadata) and a metadata with fastq file of reads that weren't trimmed
    (inserts, untrimmed_meta) = extract_inserts(joinChannel)

    // plot histogram of insert lengths
    plot_insert_histogram(inserts)

    //extract cut sites from constructs
    sites_fasta = extract_sites(joinChannel)

    input_ch
        .join(barcodes)
        .join(inserts)
        .join(sites_fasta)
        .map { meta, reads, construct, barcodes, inserts, sites ->
            // Return the meta data, barcodes, inserts, and cut sites
            [meta, barcodes, inserts, sites]
        } | set {dataChannel}

   empty_insert_histogram(dataChannel)
    // Join untrimmed_meta with input_ch, yielding a channel that contains the metadata,
    // reads, construct, flanking sequence and the fastq file of reads that weren't trimmed.
     joinChannel.join(untrimmed_meta)
        .map { meta, reads, construct, flanking, untrimmed_fq ->
            // Return the meta data, flanking sequence, and the fastq file of reads that weren't trimmed
            [meta, flanking, untrimmed_fq]
        } | set {untrimmed_with_reads}
    // extract any inserts with truncated flanking sequences
    truncated_data = extract_inserts_with_truncated_flanks(untrimmed_with_reads)

    // Generate a .tsv of the truncated insert data
    get_truncated_inserts_as_tsv(truncated_data)

    // We need the insert fasta file for the full inserts and the truncated inserts to compare them.
    data_for_plot = truncated_data.join(inserts)
         .map { meta, truncated_insert_fasta,  untrimmed_from_truncated_fq, full_inserts ->
             // Return the meta data, full inserts, truncated inserts, and untrimmed fastq file
             [meta, full_inserts, truncated_insert_fasta, untrimmed_from_truncated_fq]
         } // We also need the original sequence file to generate a list of sequence ids/lengths.
         .join(read_ch)

    // Plot the comparisons between full and truncated inserts
   plot_comparison_of_full_to_truncated_inserts(data_for_plot)

    // combine for read stats
    combined_data = input_ch
        .join(inserts)
        .join(barcodes)


    // The seq stats map is a file that maps the sample to the path of the seq stats results file for each sample.
    seq_stats_results = seq_stats(combined_data).collectFile(){
        id, tsv ->
        ["seq_stats_map.tsv", "${id}\t${params.path_prefix}${tsv}\n"]
        }

    // The sample map is a file that maps the sample to the path of the barcode count file for each sample.
    barcode_count_sample_map = barcode_counts(barcodes).collectFile(){ id, tsv ->
        ["sample_map.tsv", "${id}\t${params.path_prefix}${tsv}\n"]
    }

    // The barcode map is a file that maps the sample to the path of the barcode file for each sample.
    barcode_map =  get_barcodes_as_tsv(barcodes).collectFile(){ meta, fasta ->
        ["barcode_map.tsv", "${meta.id}\t${params.path_prefix}${fasta}\n"]
    }

    // The insert map is a file that maps the sample to the path of the insert file for each sample.
    insert_map =  get_inserts_as_tsv(inserts).collectFile(){ meta, fasta ->
        ["insert_map.tsv", "${meta.id}\t${params.path_prefix}${fasta}\n"]
    }



    // map inserts and add the dynamically generate the path to the genomic.fna file, adding it to the channel
    mapped = map_inserts(inserts | map {
	meta, seq_path -> [meta, seq_path, file(params.genomes + meta.genome + "/*_genomic.fna")]
	})

    // Add the files needed for insert coverage by dynamically creating paths based on genome.
    mapped_with_references = mapped | map {
	meta, bam, bai, stats -> [
		meta, bam]
	}

     // Here we generate various summary files and plots for all the sequences processed
     summarize_barcode_counts(barcode_count_sample_map)
     summarize_inserts(insert_map)
     summarize_barcodes(barcode_map)
     seq_summary_results = generate_seq_summary(seq_stats_results, barcode_map, vector_map, insert_map)


    // If we want to, map all reads to the donor genome
    // Add the genome .fna file to the input ch

    if (params.map_genome) {
        genome_map = map_genome(input_ch | map {
            meta, reads, construct -> [meta, reads, construct, file(params.genomes + meta.genome + "/*_genomic.fna")]
        }).collectFile(){
            meta, bam, bai, stats ->
            ["mapped_genome_map.tsv", "${meta.id}\t${params.path_prefix}${stats}\n"]
        }
        seq_stats_results.view()
        summarize_genome_mapping(genome_map, seq_stats_results)
    }
     // Plot coverage depth
     plot_depth(mapped)

}

workflow {

    log.info """
▗▄▄▖▗▄▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▄▖▗▄▄▖     ▗▄▄▖▗▄▄▄▖▗▄▄▖ ▗▄▄▄▖▗▖   ▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖ ▗▄▄▖
▐▌ ▐▌ █  ▐▌ ▐▌▐▛▚▖▐▌▐▌   ▐▌   ▐▌ ▐▌    ▐▌ ▐▌ █  ▐▌ ▐▌▐▌   ▐▌     █  ▐▛▚▖▐▌▐▌   ▐▌
▐▛▀▘  █  ▐▌ ▐▌▐▌ ▝▜▌▐▛▀▀▘▐▛▀▀▘▐▛▀▚▖    ▐▛▀▘  █  ▐▛▀▘ ▐▛▀▀▘▐▌     █  ▐▌ ▝▜▌▐▛▀▀▘ ▝▀▚▖
▐▌  ▗▄█▄▖▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▙▄▄▖▐▌ ▐▌    ▐▌  ▗▄█▄▖▐▌   ▐▙▄▄▖▐▙▄▄▖▗▄█▄▖▐▌  ▐▌▐▙▄▄▖▗▄▄▞▘

Long Read Processing and QC Pipeline
"""

 // Show help message
    print("Long read qc pipeline started and results are stored in :" + params.outdir)
    if (params.help) {
        helpMessage()
        exit 0
    }

    // If processing genome tags, the genome is not validated
    if (params.preprocess_genome_tags) {
        log.info "Preprocessing genome tags..."
        preprocess_genome_tags()
    } else{
        long_read_qc()
     }
}
