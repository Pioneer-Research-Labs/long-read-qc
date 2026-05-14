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
- The pipeline has three entry workflows implemented in `main.nf`:
  - `long_read_qc` *(default)*: the primary QC and analysis workflow that extracts barcodes and inserts, maps inserts to vector/genome, and generates per-sample and aggregate summaries.
  - `preprocess_genome_tags`: processes Tesseract multiplexed sample sheets and generates an aggregated sample sheet suitable for the main pipeline.
  - `catalog_qc`: for libraries where sequences come from a synthesized or otherwise pre-defined catalog. Maps reads to the catalog sequences instead of a reference genome. See [Catalog QC mode](#catalog-qc-mode-synthesized-libraries) below.

## Preparing construct files
- The pipeline uses SnapGene `.dna` annotations to determine flanking sequences for barcode and insert extraction. Annotate plasmids with the following feature names:
  - `BARCODEUP`, `BARCODEDN` — 5' and 3' flanking regions for the barcode
  - `INSERTUP`, `INSERTDN` — 5' and 3' flanking regions for the insert
  - `Secondary_Barcode_for_Donor_gDNA` — optional: for mixed genome libraries to denote donor-specific barcode
- Save and upload `.dna` files to your constructs S3 path ( `s3://pioneer-sequencing/constructs/`).
- Below is a schematic showing the layout of the barcode and insert features in a vector.
<img width="801" height="258" alt="Screenshot 2026-02-17 at 2 57 49 PM" src="https://github.com/user-attachments/assets/4e4fb731-90bf-4a5e-90c0-4e8772899ee6" />

- For mixed genome libraries, the `Secondary_Barcode_for_Donor_gDNA` feature can be used to denote a donor-specific barcode region that is expected to be present in reads originating from the donor genome.
This allows the pipeline to classify reads from specific genomes based on the presence of the secondary barcode. The genome-specific tags are defined in the Tesseract_Barcode_Tracking_Grid.csv and 
the pipeline searches for these tags in the reads to assign them to the correct genome. New genome tags can be added to the Tesseract_Barcode_Tracking_Grid.csv and the corresponding feature can be added to the construct .dna files to support additional genomes as needed.

- By default, the adapter lengths used by cutadapt for searching the secondary barcode is 35bp, but you can specify custom lengths with `--genome_flank_size` if needed. The
image below illustrates the placement of the adapters(in purple) for the secondary barcode in relation to the main barcode and insert features.

<img width="1082" height="272" alt="Screenshot 2026-02-26 at 12 58 32 PM" src="https://github.com/user-attachments/assets/a877a635-d4f6-4bb0-a87d-4507e9b1cd69" />


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
There are different [command line options](https://www.nextflow.io/docs/edge/cli.html) for running Nextflow.
Here we show some common ones for running the pipeline from the Insight server. Use `--help` to see the full list of options and parameters.


Log on to the Insight server, pull the latest version:

- Pull the latest pipeline (optional):

```bash
nextflow pull Pioneer-Research-Labs/long-read-qc
```

### Generating samplesheets

Use `bin/make_samplesheet.py` to generate a samplesheet from an S3 path. The script finds the FASTQ file automatically, infers the sample ID from the directory name, and uploads the samplesheet to the same S3 directory.

**Standard mode** (one sample, one reference genome):

```bash
python bin/make_samplesheet.py \
    --s3-path s3://pioneer-in-house-sequencing/ONT_minion/MY-001/ \
    --construct c.00432.dna \
    --genome H_elongata
```

**Tesseract mode** (multiplexed, genome-tagged library):

```bash
python bin/make_samplesheet.py \
    --s3-path s3://pioneer-in-house-sequencing/ONT_minion/MY-001/ \
    --construct c.00432.dna \
    --mode tesseract
```

Both commands print the S3 path of the generated samplesheet and the exact `nextflow run` command to use. Use `--id MY_SAMPLE_ID` to override the inferred sample name.

For **catalog QC mode** (synthesized library), see [Catalog QC mode](#catalog-qc-mode-synthesized-libraries) below which has its own helper script `bin/make_catalog_samplesheet.py`.

### Running the pipeline

- Run locally on the Insight server with a sample sheet entitled `samples.csv` located in the same directory:

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet samples.csv
```
- Run without Docker using your own environment (ensure all required tools are installed and in PATH):

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet samples.csv -profile local
````

- Run on AWS Batch on the Insight server with a sample sheet entitled `samples.csv` located in the same directory:

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet samples.csv -profile awsbatch
```

- Preprocess multiplexed (Tesseract) sample sheet to generate an aggregated sample sheet and run on the Insight server:

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --tesseract_samplesheet my_multiplexed_samples.csv --preprocess_genome_tags true
```

- Preprocess multiplexed (Tesseract) sample sheet on AWS Batch from the Insight server:

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --tesseract_samplesheet my_multiplexed_samples.csv --preprocess_genome_tags true -profile awsbatch
```

If you wish to clone the repo locally from your machine, you can do so and run the pipeline from your local environment.
Ensure you have access to the S3 paths for inputs and outputs and set up AWS credentials if running on AWS Batch, see these [instructions](https://www.notion.so/pioneer-labs/Running-developing-Nextflow-pipelines-locally-287c53344751809399a3f54bfc6cc337
) for configuring your local environment to run the pipeline.


Splitting large FASTQ files
- For very large FASTQ files we recommend splitting into chunked files and using the `awsbatch` profile for scalable processing.
- A helper script `bin/chunk_fastq.py` produces chunked files and an updated chunked sample sheet:

```bash
python bin/chunk_fastq.py <path to samplesheet>
```

- Run catalog QC mode (synthesized library mapped against a known sequence catalog):

```bash
nextflow run Pioneer-Research-Labs/long-read-qc --catalog_qc --catalog_samplesheet catalog_samplesheet.csv -profile awsbatch
```

---

## Catalog QC mode (synthesized libraries)

Use this mode when your FASTQ reads come from a library of synthesized sequences and you have a CSV of all designed sequences to compare against. Instead of aligning inserts to a reference genome, the pipeline converts your catalog CSV to a FASTA and maps extracted inserts against it.

All standard QC steps are preserved (quality report, vector mapping, barcode extraction, insert extraction, seq stats), and the genome alignment step is replaced with a catalog alignment.

### Catalog samplesheet format

| column | description |
|---|---|
| `id` | unique sample identifier |
| `construct` | construct `.dna` filename (under `s3://pioneer-sequencing/constructs/`) |
| `catalog` | S3 path to the synthesized sequence catalog CSV |
| `file` | S3 path to the FASTQ reads |

The catalog CSV must have these columns: `Name`, `Insert Length`, `Insert Sequence`.

### Step 1 — Generate the samplesheet

Use the helper script `bin/make_catalog_samplesheet.py` to auto-detect the FASTQ file in an S3 directory and write a ready-to-use samplesheet back to that same directory:

```bash
python bin/make_catalog_samplesheet.py \
    --s3-path s3://pioneer-in-house-sequencing/ONT_minion/HEM-003/ \
    --construct c.00631.dna \
    --catalog s3://pioneer-data/HEM-001_SUBMITTED_SEQUENCES_to_TWIST.csv
```

The script will print the S3 path of the generated samplesheet and the exact `nextflow run` command to copy-paste.

Optional flag: `--id MY_SAMPLE_ID` to override the sample ID inferred from the S3 path.

### Step 2 — Run the pipeline

```bash
nextflow run Pioneer-Research-Labs/long-read-qc \
    --catalog_qc \
    --catalog_samplesheet s3://pioneer-in-house-sequencing/ONT_minion/HEM-003/catalog_samplesheet.csv \
    -profile awsbatch
```

### Catalog QC outputs (per sample)

| file | location | description |
|---|---|---|
| `fastplong.html` | `raw_qc/` | raw read QC report |
| `mapped_vector.bam` | `primary_data/` | reads mapped to construct |
| `barcodes.fasta` | `primary_data/` | extracted barcodes |
| `inserts.fasta` | `primary_data/` | extracted inserts |
| `mapped_inserts.bam` | `primary_data/` | inserts mapped to catalog sequences |
| `mapped_inserts_stats.tsv` | `primary_data/` | flagstat summary for catalog mapping |
| `barcodes.tsv` | `primary_data/` | barcodes as TSV |
| `inserts.tsv` | `primary_data/` | inserts as TSV |
| `barcode_counts.tsv` | `summary_and_plots/` | barcode frequency counts |
| `seq_stats.tsv` | `summary_and_plots/` | read/insert/barcode counts and lengths |
| `catalog_coverage_summary.tsv` | `summary_and_plots/` | per-sequence read counts and presence/absence across the catalog |

The `catalog_coverage_summary.tsv` has columns `Name`, `Insert_Length`, `Reads_Mapped`, `Present` and reports how many reads mapped to each designed sequence.

---

### Parameters overview
- `--samplesheet <file>`: path to the sample sheet (default: `samplesheet.csv`)
- `--tesseract_samplesheet <file>`: path to Tesseract multiplexed sample sheet for preprocessing
- `--tech <str>`: sequencing tech hint (e.g., `map-ont`, `map-pb`, `map-hifi`)
- `--map_genome <bool>`: map reads to donor genome (default: `false`)
- `--preprocess_genome_tags <bool>`: run `preprocess_genome_tags` workflow instead of the main `long_read_qc` workflow
- `--error_rate <float>`: barcode search error rate (default: `0.1`)
- `--min_overlap <int>`: minimum overlap for barcode searching (default: `3`)
- `--min_bc_len` / `--max_bc_len`: minimum/maximum barcode lengths (default: `20` / `60`)
- `--cores <int>`: number of CPU cores to allocate
- `--catalog_qc <bool>`: run `catalog_qc` workflow (synthesized library mode, no reference genome; default: `false`)
- `--catalog_samplesheet <file>`: path to catalog-mode sample sheet (default: `catalog_samplesheet.csv`)
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

