# Outputs

This pipeline produces the following outputs grouped by sample id (e.g. `b.74_day0` and `s.D.1515_day0`, with summary files and plots in a `summary_and_plots` folder:
```
|-- b.74_day0
|   |-- map_genome_analysis # only if `--map_genome true`)
|   |   |-- mapped_genome.bam - bam output from mapping all reads to the genome using minimap2
|   |   |-- mapped_genome.bam.bai
|   |   `-- mapped_genome_stats.tsv - tsv output from samtools flagstats using the mapped_genome.bam file
|   |-- primary_data
|   |   |-- barcode_counts.tsv - tsv file of counts of each barcode found in the sample
|   |   |-- barcodes.fasta -fasta file of all reads with a min length of 0 in which a barcode UP/DOWN match was found
|   |   |-- barcodes.tsv - tsv file of all fastq reads with a min length of 0 in which a barcode UP/DOWN match was found
|   |   |-- flanking.gb - Genbank file of the flanking regions used for searching for barcodes, inserts, sites, and genome tags
|   |   |-- inserts.fasta - fasta file of all reads with a min length of 0 in which a insert UP/DOWN match was found
|   |   |-- inserts.tsv - tsv file of all fastq reads with a min length of 0 in which a insert UP/DOWN match was found
|   |   |-- inserts_from_truncated_flanks.fasta
|   |   |-- mapped_insert_stats.tsv - tsv output from samtools flagstats using the mapped_inserts.bam file
|   |   |-- mapped_inserts.bam - bam output from mapping inserts to the genome using minimap2
|   |   |-- mapped_inserts.bam.bai
|   |   |-- mapped_vector.bam - bam output from mapping inserts to the plasmid using minimap2
|   |   |-- mapped_vector.bam.bai
|   |   |-- mapped_vector_stats.tsv -tsv output from samtools flagstats using the mapped_vector.bam file
|   |   |-- sites.fasta -- fasta file of all reads with a min length of 0 in which a site SITEUP/DN match was found.
|   |   `-- untrimmed.fastq - fastq file of untrimmed reads in which no insert flanking region match was found.
|   |-- raw_qc
|   |   |-- fastplong.fq - qc of fastq reads
|   |   |-- fastplong.html - qc of fastq reads
|   |   `-- untrimmed_from_truncated_flanks.fastq - fastq file of untrimmed reads in which no incomplete insert flanking region match was found
|   |-- summary_and_plots
|   |   |-- counts_of_inserts_barcodes_flanking_site_sequences.csv
|   |   |-- coverage_plot.png - plot of insert length across genome based on the samtools depth output
|   |   |-- histogram_of_sites_with_barcode_no_insert.png - histogram of counts of reads in which a barcode was found but no insert
|   |   |-- insert_length_histogram.png - histogram of insert lengths for reads
|   |   |-- seq_stats.tsv - a tsv file summarizing the reads in barcode, insert, and sample reads - including num_seqs, total length, min_len, avg_len, and max_len.
|   |   `-- sites.tsv -  tsv file of all fastq reads with a min length of 0 in which a SITE UP/DOWN match was found
|   `-- truncated_analysis # only if `--analyze_truncated_flanks true`
|       |-- intact_vs_truncated_insert_lengths.png - histogram comparing insert lengths for reads in which an insert flanking region match was found vs. insert lengths for reads in which an incomplete insert flanking region match was found
|       |-- raw_counts_intact_vs_truncated_inserts.png - barplot comparing counts of reads in which an insert flanking region match was found vs. counts of reads in which an incomplete insert flanking region match was found
|       |-- read_lengths_for_trimmed_seqs.png - histogram comparing read lengths for reads in which an insert flanking region match was found vs. read lengths for reads in which an incomplete insert flanking region match was found
|       |-- truncated_inserts.tsv - tsv file of reads in which an incomplete insert flanking region match was found
|       `-- untrimmed_read_lengths.png - histogram of read lengths for untrimmed reads where an incomplete insert flanking region match was found
|-- s.D.1515_day0
`-- summary_and_plots
    |-- barcode_copy_number.png - barplot of counts of barcodes found in each sample
    |-- genome_mapping_summary.csv - csv summarizing raw_reads, reads mapped to genome, read with inserts, and reads with barcodes for each sample (if `--map_genome true`)
    |-- insert_length_barplot.png - barplot of counts of insert lengths for each sample
    `-- seq_summary.csv - csv summarizing count of raws reads, count of reads that mapped to the plasmid, the count of reads with barcodes, the count of reads with inserts, the count of reads with barcodes and inserts and the count of reads with sites for each sample. Note that the count of reads mapped to plasmid is derived from mapped_vector_stats.tsv.
```
For samtools flagstats outputs, see the [samtools documentation](https://www.htslib.org/doc/samtools.html#flagstat) for details on the fields in the output tsv files. 