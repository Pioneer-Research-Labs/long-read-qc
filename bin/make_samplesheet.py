#!/usr/bin/env python3
"""
Generate a samplesheet for the long-read-qc pipeline and upload it to the
same S3 directory as the FASTQ file.

Standard mode (maps inserts to a reference genome):
    make_samplesheet.py \\
        --s3-path s3://pioneer-in-house-sequencing/ONT_minion/MY-001/ \\
        --construct c.00432.dna \\
        --genome H_elongata \\
        [--id MY-001]

Tesseract mode (multiplexed, genome-tagged library):
    make_samplesheet.py \\
        --s3-path s3://pioneer-in-house-sequencing/ONT_minion/MY-001/ \\
        --construct c.00432.dna \\
        --mode tesseract \\
        [--id MY-001]

The script will:
  1. Find the FASTQ file at the given S3 path
  2. Infer the sample ID from the S3 path (or use --id)
  3. Write the samplesheet CSV and upload it to the same S3 directory

Then run the pipeline:
  Standard:   nextflow run Pioneer-Research-Labs/long-read-qc --samplesheet <path>
  Tesseract:  nextflow run Pioneer-Research-Labs/long-read-qc --tesseract_samplesheet <path> --preprocess_genome_tags
"""

import argparse
import csv
import io
import subprocess
import sys


FASTQ_EXTENSIONS = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')


def s3_ls(s3_dir):
    """List files directly under an S3 prefix. Returns full s3:// URIs."""
    s3_dir = s3_dir.rstrip('/') + '/'
    result = subprocess.run(
        ['aws', 's3', 'ls', s3_dir],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return []
    uris = []
    for line in result.stdout.splitlines():
        parts = line.split()
        if parts:
            uris.append(s3_dir + parts[-1])
    return uris


def find_fastq(s3_path):
    """Return the FASTQ S3 URI. Accepts a direct file path or a directory."""
    if any(s3_path.endswith(ext) for ext in FASTQ_EXTENSIONS):
        return s3_path
    files = s3_ls(s3_path)
    fastqs = [f for f in files if any(f.endswith(ext) for ext in FASTQ_EXTENSIONS)]
    if not fastqs:
        return None
    if len(fastqs) > 1:
        print(f"WARNING: Multiple FASTQ files found, using: {fastqs[0]}", file=sys.stderr)
        for f in fastqs:
            print(f"  {f}", file=sys.stderr)
    return fastqs[0]


def infer_sample_id(s3_path):
    """Derive a sample ID from the S3 path directory name.

    s3://bucket/ONT_minion/MY-001/          → 'MY-001'
    s3://bucket/ONT_minion/MY-001/reads.fq  → 'MY-001'
    """
    parts = s3_path.rstrip('/').split('/')
    for part in reversed(parts):
        if part and not any(part.endswith(ext) for ext in FASTQ_EXTENSIONS):
            return part
    return 'sample'


def main():
    parser = argparse.ArgumentParser(
        description='Generate a long-read-qc samplesheet and upload it to S3'
    )
    parser.add_argument(
        '--s3-path', required=True,
        help='S3 path to the directory containing the FASTQ, or direct path to the FASTQ file'
    )
    parser.add_argument(
        '--construct', required=True,
        help='Construct filename (e.g. c.00432.dna)'
    )
    parser.add_argument(
        '--mode', choices=['standard', 'tesseract'], default='standard',
        help='Workflow mode: standard (default, maps to reference genome) or tesseract (multiplexed, genome-tagged)'
    )
    parser.add_argument(
        '--genome',
        help='Genome ID for standard mode (e.g. H_elongata). Required when --mode standard.'
    )
    parser.add_argument(
        '--id', dest='sample_id',
        help='Sample ID (inferred from the S3 path if not provided)'
    )
    args = parser.parse_args()

    if args.mode == 'standard' and not args.genome:
        parser.error('--genome is required when --mode standard')

    # Find the FASTQ file
    fastq_path = find_fastq(args.s3_path)
    if not fastq_path:
        print(f"ERROR: No FASTQ file found at {args.s3_path}", file=sys.stderr)
        sys.exit(1)

    # Determine sample ID
    sample_id = args.sample_id or infer_sample_id(args.s3_path)

    # Output samplesheet goes in the same S3 directory as the FASTQ
    fastq_dir = '/'.join(fastq_path.split('/')[:-1])

    # Build CSV content and set output filename / run command based on mode
    buf = io.StringIO()
    writer = csv.writer(buf)

    if args.mode == 'standard':
        output_s3 = f"{fastq_dir}/samplesheet.csv"
        writer.writerow(['id', 'genome', 'construct', 'file'])
        writer.writerow([sample_id, args.genome, args.construct, fastq_path])
        run_cmd = (
            f"nextflow run Pioneer-Research-Labs/long-read-qc \\\n"
            f"    --samplesheet {output_s3}"
        )
    else:  # tesseract
        output_s3 = f"{fastq_dir}/tesseract_samplesheet.csv"
        writer.writerow(['id', 'construct', 'file'])
        writer.writerow([sample_id, args.construct, fastq_path])
        run_cmd = (
            f"nextflow run Pioneer-Research-Labs/long-read-qc \\\n"
            f"    --tesseract_samplesheet {output_s3} \\\n"
            f"    --preprocess_genome_tags"
        )

    csv_content = buf.getvalue()

    # Upload to S3 via stdin
    proc = subprocess.run(
        ['aws', 's3', 'cp', '-', output_s3],
        input=csv_content, text=True, capture_output=True
    )
    if proc.returncode != 0:
        print(f"ERROR: Failed to upload samplesheet:\n{proc.stderr}", file=sys.stderr)
        sys.exit(1)

    print(f"Samplesheet written to: {output_s3}")
    print(f"\nContents:\n{csv_content}")
    print(f"Run the pipeline with:\n  {run_cmd}")


if __name__ == '__main__':
    main()
