#!/usr/bin/env python3
"""
Generate a catalog-mode samplesheet for the long-read-qc pipeline and upload it
to the same S3 directory as the FASTQ file.

Usage:
    make_catalog_samplesheet.py \\
        --s3-path s3://pioneer-in-house-sequencing/ONT_minion/HEM-003/ \\
        --construct c.00631.dna \\
        --catalog s3://pioneer-data/HEM-001_SUBMITTED_SEQUENCES_to_TWIST.csv \\
        [--id HEM-003]

The script will:
  1. Find the FASTQ file at the given S3 path (or use it directly if a file path is given)
  2. Infer the sample ID from the S3 path (or use --id if provided)
  3. Write a catalog_samplesheet.csv and upload it to the same S3 directory

Then run the pipeline with:
    nextflow run Pioneer-Research-Labs/long-read-qc \\
        --catalog_qc \\
        --catalog_samplesheet <printed S3 path>
"""

import argparse
import csv
import io
import subprocess
import sys


FASTQ_EXTENSIONS = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')


def s3_ls(s3_dir):
    """List filenames directly under an S3 prefix. Returns full s3:// URIs."""
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
    """Derive a sample ID from the S3 path.

    For a directory path like  s3://bucket/ONT_minion/HEM-003/  → 'HEM-003'
    For a file path like       s3://bucket/ONT_minion/HEM-003/combined.fastq → 'HEM-003'
    """
    parts = s3_path.rstrip('/').split('/')
    # Walk backwards; skip the filename if it looks like a fastq
    for part in reversed(parts):
        if part and not any(part.endswith(ext) for ext in FASTQ_EXTENSIONS):
            return part
    return 'sample'


def main():
    parser = argparse.ArgumentParser(
        description='Generate a catalog-mode samplesheet and upload it to S3'
    )
    parser.add_argument(
        '--s3-path', required=True,
        help='S3 path to the directory containing the FASTQ, or direct path to the FASTQ file'
    )
    parser.add_argument(
        '--construct', required=True,
        help='Construct filename (e.g. c.00631.dna)'
    )
    parser.add_argument(
        '--catalog', required=True,
        help='S3 path to the synthesized sequence catalog CSV'
    )
    parser.add_argument(
        '--id', dest='sample_id',
        help='Sample ID (inferred from the S3 path if not provided)'
    )
    args = parser.parse_args()

    # Find the FASTQ file
    fastq_path = find_fastq(args.s3_path)
    if not fastq_path:
        print(f"ERROR: No FASTQ file found at {args.s3_path}", file=sys.stderr)
        sys.exit(1)

    # Determine sample ID
    sample_id = args.sample_id or infer_sample_id(args.s3_path)

    # Output samplesheet goes in the same S3 directory as the FASTQ
    fastq_dir = '/'.join(fastq_path.split('/')[:-1])
    output_s3 = f"{fastq_dir}/catalog_samplesheet.csv"

    # Build CSV content
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(['id', 'construct', 'catalog', 'file'])
    writer.writerow([sample_id, args.construct, args.catalog, fastq_path])
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
    print(f"Run the pipeline with:")
    print(f"  nextflow run Pioneer-Research-Labs/long-read-qc \\")
    print(f"    --catalog_qc \\")
    print(f"    --catalog_samplesheet {output_s3}")


if __name__ == '__main__':
    main()
