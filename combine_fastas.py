#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

def combine_fasta(files, output):
    seen_ids = set()
    unique_records = []

    for file in files:
        with open(file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id not in seen_ids:
                    seen_ids.add(record.id)
                    unique_records.append(record)

    with open(output, "w") as out_handle:
        SeqIO.write(unique_records, out_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description='Combine multiple FASTA files into one, removing redundant entries.')
    parser.add_argument('--files', type=str, required=True, nargs='+', help='Paths to the input FASTA files')
    parser.add_argument('--output', type=str, required=True, help='Path to the output FASTA file')
    args = parser.parse_args()

    combine_fasta(args.files, args.output)

if __name__ == "__main__":
    main()