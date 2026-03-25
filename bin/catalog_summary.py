#!/usr/bin/env python3
"""
Generate a per-sequence coverage summary for a known-sequence catalog.

Joins the catalog CSV (Name, Insert Length, Insert Sequence) with
samtools idxstats output to report how many reads mapped to each designed sequence.
"""

import csv
import sys


def main():
    catalog_csv = sys.argv[1]
    idxstats_file = sys.argv[2]
    output_file = sys.argv[3]

    # Parse samtools idxstats: name, seq_length, mapped_reads, unmapped_reads
    reads_per_seq = {}
    with open(idxstats_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[0] != '*':
                reads_per_seq[parts[0]] = int(parts[2])

    # Read catalog CSV and join
    with open(catalog_csv) as f:
        sample = f.read(2048)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters='\t,')
        reader = csv.DictReader(f, dialect=dialect)
        rows = list(reader)

    total = len(rows)
    found = 0

    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['Name', 'Insert_Length', 'Reads_Mapped', 'Present'])
        for row in rows:
            name = row['Name'].strip()
            length = row['Insert Length'].strip()
            mapped = reads_per_seq.get(name, 0)
            present = 'yes' if mapped > 0 else 'no'
            writer.writerow([name, length, mapped, present])
            if mapped > 0:
                found += 1

    pct = f"{found / total * 100:.1f}" if total > 0 else "N/A"
    print(f"Catalog coverage: {found}/{total} sequences detected ({pct}%)")


if __name__ == '__main__':
    main()
