#!/usr/bin/env python3
"""Convert a synthesized library CSV (Name, Insert Length, Insert Sequence) to FASTA."""

import csv
import sys


def main():
    csv_file = sys.argv[1]
    fasta_file = sys.argv[2]

    with open(csv_file) as f:
        sample = f.read(2048)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters='\t,')
        reader = csv.DictReader(f, dialect=dialect)
        with open(fasta_file, 'w') as out:
            for row in reader:
                name = row['Name'].strip()
                seq = row['Insert Sequence'].strip()
                out.write(f">{name}\n{seq}\n")


if __name__ == '__main__':
    main()
