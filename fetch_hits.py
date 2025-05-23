#!/bin/python
import sys
import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='Fetches fasta hits from database file.')
    parser.add_argument('--database', help='Path to individual GTDB fastas, named with species ID', required=True)
    parser.add_argument('--hits', help='Hits file', required=True)
    args = parser.parse_args()

    database = args.database
    hits_file = args.hits
    output_file = os.path.basename(hits_file).rstrip('.txt')+'.fa'
    imported_proteins = {}
    gene_count = 0
    hits_count = 0    
    fastas = os.listdir(database)
    total = len(fastas)
    i = 0
    for fasta in fastas:
        i+=1
        update_counter(i, total)
        fasta = os.path.join(database, fasta)
        imported_proteins, gene_count = import_fasta(fasta, imported_proteins, gene_count)
    output(hits_file, output_file, imported_proteins, hits_count, gene_count)
    
def update_counter(count, total):
    sys.stdout.write("\rProgress: {}/{} fastas loaded".format(count, total))
    sys.stdout.flush()
    
def import_fasta(fasta, imported_proteins, gene_count):
    with open(fasta, 'r') as fasta_file:  
        sequence = None
        for line in fasta_file:
            if '|' in line:
                continue
            if '>' in line:
                gene_count += 1
                if sequence is not None:
                    imported_proteins[protein] = [species, sequence]
                    sequence = None
                protein = (line.split(' ', 1)[0].lstrip('>'))
                species = os.path.splitext(os.path.basename(fasta))[0].rstrip('_protein')
            elif (len(line) != 0) and (sequence is not None):
                sequence += (line.rstrip('\n'))
            elif (len(line) != 0) and (sequence is None):
                sequence = (line.rstrip('\n'))
            imported_proteins[protein] = [species, sequence]
    return imported_proteins, gene_count
            
def output(hits_file, output_file, imported_proteins, hits_count, gene_count):
    with open(hits_file, 'r') as hits_file:
        with open(output_file, 'w') as output_file:
            for line in hits_file:
                hits_count += 1
                protein = line.rstrip('\n')
                species = imported_proteins[protein][0]
                sequence = imported_proteins[protein][1]
                output_file.write('>'+protein+' ['+species+']\n')
                output_file.write(sequence+'\n')
    print('\n'+str(hits_count)+' hits out of '+str(gene_count)+' total')

if __name__ == '__main__':
    main()