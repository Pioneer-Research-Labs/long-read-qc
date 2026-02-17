# long-read-qc

Long-read processing and QC pipeline for Pioneer Research Labs.

This repository contains a Nextflow-based pipeline that extracts barcodes and inserts from long-read sequencing data, maps inserts to vector and genome references, and generates QC reports and summary tables/plots.

## Key features
- Modular Nextflow workflows for extraction, mapping, QC and summarization.
- Supports preprocessing of genome tags and truncated flank analysis.
- Produces per-sample results and aggregate summary maps for downstream reporting.

## Repository layout
- `main.nf` — pipeline entrypoint and workflow orchestration.
- `modules/` — Nextflow modules implementing atomic tasks (barcode extraction, mapping, plotting, etc.).
- `bin/` — helper scripts used by modules and workflow (aggregation, conversions, etc.).
- `assets/` — validation schemas, templates and example resources.
- `tests/` — Nextflow tests for helper scripts and modules.
- `test_data/` — small example inputs and expected artifacts for local testing.

## Requirements
- Java 11+ (for Nextflow runtime)
- Nextflow (https://www.nextflow.io)
- Docker or Singularity (optional; pipeline can run locally or inside containers)
- Standard bioinformatics tools used inside modules (see `modules/*` for per-module requirements). Many modules expect standard command-line tools (samtools, minimap2, cutadapt, seqkit, bedtools, etc.).

## Quick install
1. Install Java and Nextflow per official docs.
2. (Optional) Install Docker or Singularity if you plan to run with containers.
3. Clone this repository and change into it:

```bash
git clone https://github.com/Pioneer-Research-Labs/long-read-qc.git
cd long-read-qc
```

## High-level workflow
- The pipeline has two main entry workflows implemented in `main.nf`:
  - `preprocess_genome_tags`: processes Tesseract multiplexed sample sheets and generates an aggregated sample sheet suitable for the main pipeline.
  - `long_read_qc`: the primary QC and analysis workflow that extracts barcodes and inserts, maps inserts to vector/genome, and generates per-sample and aggregate summaries.

## Preparing construct files
- The pipeline uses SnapGene `.dna` annotations to determine flanking sequences for barcode and insert extraction. Annotate plasmids with the following feature names:
  - `BARCODEUP`, `BARCODEDN` — 5' and 3' flanking regions for the barcode
  - `INSERTUP`, `INSERTDN` — 5' and 3' flanking regions for the insert
  - `SITEUP`, `SITEDN` — 5' and 3' flanking regions surrounding insertion site(s)
  - `Secondary_Barcode_for_Donor_gDNA` — optional: for mixed genome libraries to denote donor-specific barcode
- Save and upload `.dna` files to your constructs S3 path ( `s3://pioneer-sequencing/constructs/`).

## Sample sheet format
- A CSV describing samples is required. Minimum columns (for standard runs): `id`, `genome`, `construct`, `file`.
  - `id`: a unique sample identifier
  - `genome`: genome id (used to find reference `*_genomic.fna` under your genome path in S3 - `s3://pioneer-sequencing/pioneer-genomes/`)
  - `construct`: plasmid name (filename of `.dna` under the construct path in `s3://pioneer-sequencing/constructs/`)
  - `file`: path to the sample reads in S3 path like `s3://.../sample.fastq.gz`)
- See an example sample sheet in `documentation/example_samplesheet.csv`.
- For mixed-genome libraries the `genome` column may be omitted or left blank; the pipeline supports this mode when using the Tesseract preprocessing workflow.
- Validation schemas are included in `assets/` (e.g. `samplesheet_validation_schema.json`) and are used to validate sample sheets when the pipeline starts in Seqera.

## Running the pipeline from Seqera
- The pipeline is designed to run on Seqera's AWS Batch infrastructure, but can also be run locally with Docker or Singularity profiles. 
Follow the instructions for running on Seqera [here](https://www.notion.so/pioneer-labs/Infrastructure-149c5334475180d18f8fecfaf3eac640?source=copy_link#216c53344751805f9e82d3e28d49a7cf) or
follow the local run instructions below.


## Running the pipeline locally
When running the pipeline locally, you'll need to authenticate with AWS to allow the pipeline to read/write from S3. 
You can do this by [configuring the AWS CLI with your credentials](https://www.notion.so/pioneer-labs/Running-developing-Nextflow-pipelines-locally-287c53344751809399a3f54bfc6cc337#287c53344751803c8f11d22a4e3f234c).

```bash
Running the pipeline
- Pull the latest pipeline (optional):

```bash
nextflow pull Pioneer-Research-Labs/long-read-qc
```
 
- Run locally with Docker profile (example):

```bash
nextflow run Pioneer-Research-Labs/long-read-qc -profile standard --samplesheet samples.csv
```

- Run on AWS Batch (example):

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet samples.csv -profile awsbatch
```

- Preprocess multiplexed (Tesseract) sample sheet and generate an aggregated sample sheet:

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --tesseract_samplesheet my_multiplexed_samples.csv --preprocess_genome_tags true
```

Splitting large FASTQ files
- For very large FASTQ files we recommend splitting into chunked files and using the `awsbatch` profile for scalable processing.
- A helper script `bin/chunk_fastq.py` produces chunked files and an updated chunked sample sheet:

```bash
python bin/chunk_fastq.py <path to samplesheet>
```

### Parameters overview
- `--samplesheet <file>`: path to the sample sheet (default: `samplesheet.csv`)
- `--tesseract_samplesheet <file>`: path to Tesseract multiplexed sample sheet for preprocessing
- `--tech <str>`: sequencing tech hint (e.g., `map-ont`, `map-pb`, `map-hifi`)
- `--map_genome <bool>`: map reads to donor genome (default: `false`)
- `--preprocess_genome_tags <bool>`: run `preprocess_genome_tags` workflow instead of main
- `--error_rate <float>`: barcode search error rate (default: `0.1`)
- `--min_overlap <int>`: minimum overlap for barcode searching (default: `3`)
- `--min_bc_len` / `--max_bc_len`: minimum/maximum barcode lengths (default: `20` / `60`)
- `--cores <int>`: number of CPU cores to allocate
- Use `--help` or inspect `main.nf`'s `helpMessage()` for the full list and defaults.

### Outputs
- The output is stored in a directory in S3 (s3://pioneer-analysis/long_read_qc/prod/).
The directory name the sample sheet file name prefixed with today's date (yyyy_MM_dd_HH_mm_ss). 
For example, a sample sheet entitled `my_samples.csv` will result in an output directory `2025_03_04_14_37_19_my_samples`.
In the output directory, a folder is created for each sample containing task outputs, summary TSVs, plots, and intermediate files as applicable.
- For more details on outputs, see `OUTPUTS.md` in the documentation folder.

### Workflow diagram
- See the workflow diagram in the documentation folder (`documentation/workflow_diagram.png`) for a visual overview of the main workflow and its modules.

## Testing
- Install `nf-test`(https://www.nf-test.com/)
- Use the provided `test_data/` to run local smoke tests

```bash
nf-test test tests
```

## Developer notes
- The pipeline is modular: each module in `modules/` implements a small Nextflow process or subworkflow. Open a module to see tool-specific commands, containers, and resource hints.
- Helper scripts live in `bin/` (e.g. `aggregate_sample_sheets.py`, `chunk_fastq.py`) and can be run standalone.
- `params.path_prefix` is used to add prefixes (for S3 or other storage) to outputs in collection maps; verify this if your outputs end up under unexpected locations. This
and other parameters can be set in `nextflow.config` or overridden on the command line.

## Troubleshooting
- Missing tools: inspect the failing task logs under `work/` (task stderr/stdout) and install the required tool or run with an appropriate container profile.
- File not found / permission errors: verify sample sheet paths and S3 access. Use absolute paths or S3 URIs, and ensure IAM permissions for S3 when running on AWS.
- Resume a failed run using Nextflow's resume feature by re-running the same command with `-resume`:

```bash
nextflow run main.nf -resume --samplesheet samples.csv
```

- For difficult failures, include the Nextflow run command, Nextflow version (`nextflow -v`), Java version (`java -version`), and the relevant `work/<task>` stderr log when opening an issue.

## Testing and CI
- Continuous integration is set up with GitHub Actions to run `nf-test` on the `tests/` directory for each push and pull request. This ensures that helper scripts and modules are functioning as expected with the provided test data.
- To add new tests, create a new test file in `tests/` (e.g. `test_new_module.nf`) and add corresponding test data in `test_data/`. The test should include a minimal Nextflow script that calls the module or script being tested, along with assertions to verify expected outputs.

## Contributing
- Please open issues and PRs. Include a minimal reproducer (use `test_data/` if possible) and describe expected vs actual behavior.

## Contact
- For help or questions, provide the pipeline run command and a snippet of failing task logs (found under `work/`).

