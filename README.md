# Analysis pipeline for long reads

Possible names:
* Extractify
* Genomicode
* Extract-a-seq
* FATTSeq (Find All The Things)
* BubbleSeq

## Usage

### Prep construct file

To determine the regions of the read that have the barcode and insert the pipeline uses flanking regions for each part.  These are provided by way of annotations in a SnapGene `.dna` file.  To create these, add annotations named `BARCODEUP`, `BARCODEDN`, `INSERTUP`, and `INSERTDN` for the 5' and 3' flanking regions for the barcode and insert.  Ideally these will be between 100-30 bp long but the exact length is not particularly important here.  Generally don't go below 20 bp but that's not a fixed limit.  Once you've added these annotations with the exact names as above the save the file and upload to the Jupyter server somewhere you can find it.

### Create sample sheet

To run the pipeline a csv file describing the samples must be provided.  At minimum this file needs to have three columns titled, `id`, `genome`, `construct`, and `file`.  The `id` column should be a unique sample id that is simple and not too long. `genome` is the genome id to use for mapping inserts and identifying genes.  See below for current list.  `construct` is the path to a SnapGene plasmid map with annotations described above.  Note this can also be an S3 URI for a file located in one of Pioneer's S3 buckets.  The `file` column is the path to the raw fastq file for that sample.  This can be a local file but it's highly recommended to have the raw data in an S3 bucket and use the S3 URI here.  

Once you've created this file (probably easiest to do this on your local machine, not the server), you can upload it to the server in the working directory you'd like to use for the run.

#### Current genomes

- B_subtilis
- C_psychrerythraea
- G_obscurus
- H_elongata
- P_halocryophilus
- R_radiotolerans

### Run the pipeline

Pull the latest version of the pipeline:

```
nextflow pull Pioneer-Research-Labs/long-read-qc
```

Start the run:

```
nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet mysample.csv
```

For low-depth QC samples this should run in less than a minute

### Results

The output is stored in a directory called `results`.  This can be changed by providing the `--outdir` paramater, i.e., `--outdir my_awesome_run`.  In the output directory a folder is created for each sample id and the raw output files are stored here.  

To run the report and perform additional analysis if required the `results` directory contains a Jupyter notebook called `report.ipynb`.  On the server you can open this directly and run it to get all the sample summaries and plots.  Please note that you may want to restart the Jupyter Kernel if you've run other notebooks recently with the same name.  The report notebook will aggregate the results from the samples, write out aggregated files for future analysis, and create some summary plots and tables.  Feel free to modify and extend the analysis as requried.


### Pipeline options

#### General

`--samplesheet`: Three column csv file with columns titled "id", "construct", and "file" (default: "samplesheet.csv")

`--outdir`: Name of output directory to save the results (default: "results")

`--rotate_anchor`: vector sequence to use for rotating and orienting reads (default: rrnB terminator, "TTTATCTGTTGTTTGTCGGTGAACGCTCTC")

#### Genome files for mapping inserts

`--ref_fa`: Fasta file of reference genome (default: "$HOME/shared/genomes/H_elongata/H_elongata_contigs.fna")

`--ref_gff`: GFF file with gene annotations (default: "$HOME/shared/genomes/H_elongata/H_elongata_genes.gff")

`--ref_bed`: BED file with gene annotations (default: "$HOME/shared/genomes/H_elongata/H_elongata_genes.bed")


#### Barcode and insert searching parameters

`--error_rate`: Error rate for matching flanking regions (default: 0.1)

`--min_overlap`: Minimum overlap to identify an flanking region (default: 3)


#### Barcode filtering

`--min_bc_len`: Minimum barcode length to include in output (default: 20)

`--max_bc_len`: Maximium barcode length to include in output (default: 60)

## About

This is a fairly straightforward pipeline that uses `cutadapt` to locate barcodes and inserts on the reads.  It first uses the [`rotate`](https://github.com/richarddurbin/rotate) to orient the reads using a provided anchor sequence present on the plasmid.  This also filters out reads that don't have plasmid sequence so best to use an anchor sequence that is close to the barcode/insert regions.  

Using the provided `.dna` Snapgene file the annotations are extracted and flanking sequences saved in a format for `cutadapt`.   Then the reads are run through `cutadapt` twice, once to extract the barcode and once to extract the insert.  Cutadapt has a feature called "linked adapters"  which will uses the provided flanking sequence to locate the 5' and 3' flanking regions and then it trims off this sequence from the read and everything outside it, leaving just the barcode or insert.

After counting unique barcodes the inserts are mapped to the provided genome and overlaps with genes and the genome are assesed using a variety of bedtools commands


