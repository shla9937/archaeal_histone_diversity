#!/usr/bin/env python
import argparse
import re
import os


def main():
    p = argparse.ArgumentParser(description="Annotate FASTA headers in file B with genome names from FASTA A.")
    p.add_argument('-s', required=True, help="FASTA A (contains accession with [GENOME] in header)")
    p.add_argument('-i', required=True, help="FASTA B (headers to annotate)")
    p.add_argument('-o', default='cluster_1_genomes.fa', help="Output FASTA path")
    args = p.parse_args()

    accession_map = build_accession_map(args.s)
    count = annotate_fasta(args.i, accession_map, args.o)
    print(f"Wrote {args.o} ({count} headers processed). Found mappings for {len(accession_map)} accessions.")

def build_accession_map(fasta_path):
    m = {}
    with open(fasta_path, 'r') as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            # match accession and optional [GENOME]
            mo = re.match(r'^>(\S+)(?:\s+\[([^\]]+)\])?', line)
            if mo:
                acc = mo.group(1)
                genome = mo.group(2)
                if genome:
                    m[acc] = genome
    return m

def annotate_fasta(input_fa, accession_map, out_path):
    written = 0
    with open(input_fa, 'r') as inf, open(out_path, 'w') as outf:
        for line in inf:
            if line.startswith('>'):
                mo = re.match(r'^>(\S+)(.*)', line.rstrip('\n'))
                acc = mo.group(1)
                rest = mo.group(2).strip()
                # if already has bracketed genome info, leave as-is
                if rest and '[' in rest and ']' in rest:
                    outf.write(line)
                else:
                    if acc in accession_map:
                        outf.write(f'>{acc} [{accession_map[acc]}]\n')
                    else:
                        outf.write(line)
                written += 1
            else:
                outf.write(line)
    return written


if __name__ == '__main__':
    main()