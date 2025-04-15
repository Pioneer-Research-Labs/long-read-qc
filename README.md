# Analysis pipeline for long reads

## Usage

### Prep construct file

To determine the regions of the read that have the barcode and insert the pipeline uses flanking regions for each part.  
These are provided by way of annotations in a SnapGene `.dna` file.  To create these, add annotations named `BARCODEUP`, 
`BARCODEDN`, `INSERTUP`, and `INSERTDN` for the 5' and 3' flanking regions for the barcode and insert.  Ideally these 
will be between 100-30 bp long but the exact length is not particularly important here.  Generally don't go below 20 bp
but that's not a fixed limit.  Once you've added these annotations with the exact names as above the save the file and
upload to the S3 path at s3://pioneer-sequencing/constructs/.

### Create sample sheet

To run the pipeline a csv file describing the samples must be provided.  At minimum this file needs to have four columns 
titled, `id`, `genome`, `construct`, and `file`.  The `id` column should be a unique sample id that is simple and not too 
long. `genome` is the genome id to use for mapping inserts and identifying genes.  See below for current list of genomes.  
`construct` is the name of the plasmid. This should be the plasmid name without the path to the S3 location.  The `file` 
column is the S3 path to the raw fastq file for that sample. The only column in the sample sheet that contains a full S3
path is the file column.  Here's an example:

| id                                 | genome     | construct   | file                                                             |
|------------------------------------|------------|-------------|------------------------------------------------------------------|
| HE_Gateway_Lib1_20241217           | H_elongata | c.00323.dna | s3://pioneer-scratch/dummy/dummy_seqs/5HKJT6_3_sample_3.fastq.gz |
| Gateway_EmptyBarcodedLib2_20241217 | H_elongata | c.00323.dna | s3://pioneer-scratch/dummy/dummy_seqs/5HKJT6_2_sample_2.fastq.gz |
| HE_USER-NOgapfill_Lib1_20241220    | H_elongata | c.00215.dna | s3://pioneer-scratch/dummy/dummy_seqs/CJ39YZ_2_sample_2.fastq.gz |


Once you've created this file (probably easiest to do this on your local machine, not the server), you can upload it 
to the server in the working directory you'd like to use for the run.

#### Current genomes available in the pipeline (s3://pioneer-data/genomes/)

- B_subtilis
- C_psychrerythraea
- G_obscurus
- H_elongata
- P_halocryophilus
- R_radiotolerans

### Large sequencing files
For large sequencing files we recommend splitting the files into smaller chunks and using the `awsbatch` pipeline profile.
To split the files, check out the repo and run the following command:

```bash
python bin/chunk_fastq.py <path to samplesheet>
```
The script takes a while to run since it first determines how many reads are present, and then calculates the number of
reads per 1GB of data.  The file is split into 1GB chunks and the chunks uploaded to S3. A new sample sheet is created
with the same name as the input file but with `_chunked` appended to the end.  
This file will have the same structure as the original sample sheet but with the `file` column updated to point to the 
chunked files loaded to S3 and the `id` column updated with sample and chunk name. 

### Run the pipeline

Before running the pipeline make sure you have the latest version of the pipeline.  You can do this by running the following command:

```
nextflow pull Pioneer-Research-Labs/long-read-qc
````

Start the run locally:

```
nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet mysample.csv 
```

Start the run on AWS Batch:

```
nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet mysample.csv -profile awsbatch
```

For low-depth QC samples this should run in less than a minute. For high-depth PacBio samples this will take an hour or more.

### Results

The output is stored in a directory locally and in S3 (s3://pioneer-analysis/long_read_qc/prod/). The directory name
the sample sheet file name prefixed with today's date (yyyy_MM_dd_HH_mm_ss). For example, a sample sheet entitled `my_samples.csv` will
result in an output directory `2025_03_04_14_37_19_my_samples`. In the output directory, a folder is created for each sample id and the raw output 
files are stored here.  

In addition to outputs for each sample, the pipeline aggregates results across all samples and generates the following plots:
- `barcode_copy_number.png` - A histogram of the copy number of each unique barcode
- `barcode_length_distribution.png` - A violin plot of the length distribution of the barcodes
- `barcode_proportions.png` - A bar chart showing the proportion of true/empty/other detected barcodes
- `full_genes_per_fragment.png` - A histogram of the number of full genes per insert
- `insert_length_distribution.png` - A violin plot of the length distribution of the inserts
- `partial_genes_per_fragment.png` - A histogram of the number of partial genes per insert

A file (seq_summary.csv) summarizing the results for each sample is also generated in the output directory.

The legacy tool to generate reports and plots is still in place.  To run the report and perform additional analysis 
if required the `results` directory contains a Jupyter notebook called `report.ipynb`.  On the server you can open 
this directly and run it to get all the sample summaries and plots.  Please note
that you may want to restart the Jupyter Kernel if you've run other notebooks recently with the same name.  The report 
notebook will aggregate the results from the samples, write out aggregated files for future analysis, and create some summary 
plots and tables.  Feel free to modify and extend the analysis as requried.


### Pipeline options

```
▗▄▄▖▗▄▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▄▖▗▄▄▖     ▗▄▄▖▗▄▄▄▖▗▄▄▖ ▗▄▄▄▖▗▖   ▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖ ▗▄▄▖
▐▌ ▐▌ █  ▐▌ ▐▌▐▛▚▖▐▌▐▌   ▐▌   ▐▌ ▐▌    ▐▌ ▐▌ █  ▐▌ ▐▌▐▌   ▐▌     █  ▐▛▚▖▐▌▐▌   ▐▌   
▐▛▀▘  █  ▐▌ ▐▌▐▌ ▝▜▌▐▛▀▀▘▐▛▀▀▘▐▛▀▚▖    ▐▛▀▘  █  ▐▛▀▘ ▐▛▀▀▘▐▌     █  ▐▌ ▝▜▌▐▛▀▀▘ ▝▀▚▖
▐▌  ▗▄█▄▖▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▙▄▄▖▐▌ ▐▌    ▐▌  ▗▄█▄▖▐▌   ▐▙▄▄▖▐▙▄▄▖▗▄█▄▖▐▌  ▐▌▐▙▄▄▖▗▄▄▞▘

Long Read Processing and QC Pipeline          


Usage: nextflow run Pioneer-Research-Labs/long-read-qc -latest

Options:
--samplesheet <file>      Path to the sample sheet (default: samplesheet.csv)
--outdir <dir>            Output directory (default: results)
--tech <str>              Sequencing technology, map-ont/map-pb/map-hifi (default: map-ont)
--error_rate <float>      Error rate for barcode searching (default: 0.1)
--min_overlap <int>       Minimum overlap for barcode searching (default: 3)
--min_bc_len <int>        Minimum barcode length (default: 20)
--max_bc_len <int>        Maximum barcode length (default: 60)
--meta_ovlp <int>         Overlap bp for sourmash (default: 1000)
--sourmash_db <file>      Path to the sourmash database (default: /srv/shared/databases/sourmash/gtdb-rs220-k21.zip)
--taxonomy <file>         Path to the taxonomy database (default: /srv/shared/databases/sourmash/gtdb-rs220.lineages.sqldb)
--cores <int>             Number of cores to use (default: 4)

Profiles:
-profile standard              Run pipeline locally with Docker
-profile awsbatch              Run pipeline on AWS Batch
```

## About

This is a fairly straightforward pipeline that uses `cutadapt` to locate barcodes and inserts on the reads. 

Using the provided `.dna` Snapgene file the annotations are extracted and flanking sequences saved in a format for `cutadapt`.   Then the reads are run through `cutadapt` twice, once to extract the barcode and once to extract the insert.  Cutadapt has a feature called "linked adapters"  which will uses the provided flanking sequence to locate the 5' and 3' flanking regions and then it trims off this sequence from the read and everything outside it, leaving just the barcode or insert.

After counting unique barcodes the inserts are mapped to the provided genome and overlaps with genes and the genome are assesed using a variety of bedtools commands


## Developing

### Adding genome files
TODO
