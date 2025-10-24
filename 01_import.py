#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
import statistics
from statistics import variance
import numpy as np
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description='Import fasta of proteins, output df for downstream.')
    parser.add_argument('--fa', required=True, type=str,
                        help='FASTA with protein sequences, species ID should be included in entry.')
    args = parser.parse_args()
    imports = import_file(str(args.fa))
    att_proteins, total, excluded = attribute_proteins(imports)
    proteins = standardize_proteins(att_proteins)
    variances(proteins, total, excluded)
    write_df(proteins)
    
def import_file(input_file):
    names = []
    species = []
    sequences = []
    lengths = []
    sequence = None
    gene_count = 0   
    with open(input_file, 'r') as fastafile:
        for line in fastafile:
            if '|' in line:
                continue
            if '>' in line:
                gene_count += 1
                if sequence is not None:
                    lengths.append(len(sequence))
                    sequences.append(sequence)
                    sequence = None
                names.append(line.split(' [')[0])
                try:
                    species.append(line.split('[')[1].split(']')[0])
                except:
                    species.append('Null')
            elif (len(line) != 0) and (sequence is not None):
                sequence += (line.rstrip('\n').rstrip('*'))
            elif (len(line) != 0) and (sequence is None):
                sequence = (line.rstrip('\n').rstrip('*'))
        lengths.append(len(sequence))
        sequences.append(sequence)
        proteins = pd.DataFrame({'Species ID': species, 'Sequence': sequences, 'Length': lengths}, index = names)
    return proteins

def attribute_proteins(proteins):
    excluded_proteins = 0
    counter = 0
    total = len(proteins)
    for protein in proteins.index:
        counter += 1
        sequence = proteins.loc[protein, 'Sequence']
        update_counter(counter, total)
        if any(substr in sequence for substr in ['X', 'O', 'U', 'B', 'Z', 'partial']):
            excluded_proteins += 1
            proteins.drop(protein, axis=0, inplace=True)
            continue
        else:
            X = ProtParam.ProteinAnalysis(sequence)
            proteins.at[protein, 'Aromaticity'] = X.aromaticity()
            proteins.at[protein, 'pI'] = X.isoelectric_point()
            proteins.at[protein, 'Helix'] = X.secondary_structure_fraction()[0]
            proteins.at[protein, 'Turn'] = X.secondary_structure_fraction()[1]
            proteins.at[protein, 'Sheet'] = X.secondary_structure_fraction()[2]
            proteins.at[protein, 'Instability'] = X.instability_index()
            proteins.at[protein, 'Flexibility'] = sum(X.flexibility() / proteins.loc[protein, 'Length'])
            proteins.at[protein, 'Gravy'] = X.gravy()
    return proteins, total, excluded_proteins

def update_counter(count, total):
    sys.stdout.write("\rProgress: {}/{} proteins processed".format(count, total))
    sys.stdout.flush()

def standardize_data(data_list):
    standard_list = []
    list_stdev = statistics.stdev(data_list)
    list_mean = statistics.mean(data_list)
    for i in data_list:
        new_i = (i - list_mean) / list_stdev
        standard_list.append(new_i)
    return standard_list

def unstandardize_data(data_list, coords):
    unstandard_list = []
    list_stdev = statistics.stdev(data_list)
    list_mean = statistics.mean(data_list)
    for i in coords:
        new_i = (i * list_stdev) + list_mean
        unstandard_list.append(new_i)
    return unstandard_list

def standardize_proteins(proteins):
    proteins['New Length'] = standardize_data(proteins['Length'])
    proteins['New Aromaticity'] = standardize_data(proteins['Aromaticity'])
    proteins['New pI'] = standardize_data(proteins['pI'])
    proteins['New Helix'] = standardize_data(proteins['Helix'])
    proteins['New Turn'] = standardize_data(proteins['Turn'])
    proteins['New Sheet'] = standardize_data(proteins['Sheet'])
    proteins['New Instability'] = standardize_data(proteins['Instability'])
    proteins['New Flexibility'] = standardize_data(proteins['Flexibility'])
    proteins['New Gravy'] = standardize_data(proteins['Gravy'])
    return proteins

def variances(proteins, total, excluded):
    with open('../outputs/protein_info.txt', 'w') as f:
        print('\n')
        print('Writing the following to ../outputs/protein_info.txt:')
        f.write('Total: '+str(total) + '\n')
        print('Total: '+str(total))
        f.write('Excluded: '+str(excluded) + '\n')
        print('Excluded: '+str(excluded))
        f.write('Included: '+str(len(proteins)) + '\n')
        print('Included: '+str(len(proteins)))
        f.write('Variance length: ' + str(variance(proteins['Length'])) + '\n')
        print('Variance length: '+str(variance(proteins['Length'])))
        f.write('Variance aromaticity: ' + str(variance(proteins['Aromaticity'])) + '\n')
        print('Variance aromaticity: '+str(variance(proteins['Aromaticity'])))
        f.write('Variance pI: ' + str(variance(proteins['pI'])) + '\n')
        print('Variance pI: '+str(variance(proteins['pI'])))
        f.write('Variance helix: ' + str(variance(proteins['Helix'])) + '\n')
        print('Variance helix: '+str(variance(proteins['Helix'])))
        f.write('Variance turn: ' + str(variance(proteins['Turn'])) + '\n')
        print('Variance turn: '+str(variance(proteins['Turn'])))
        f.write('Variance sheet: ' + str(variance(proteins['Sheet'])) + '\n')
        print('Variance sheet: '+str(variance(proteins['Sheet'])))
        f.write('Variance instability: ' + str(variance(proteins['Instability'])) + '\n')
        print('Variance instability: '+str(variance(proteins['Instability'])))
        f.write('Variance flexibility: ' + str(variance(proteins['Flexibility'])) + '\n')
        print('Variance flexibility: '+str(variance(proteins['Flexibility'])))
        f.write('Variance gravy: ' + str(variance(proteins['Gravy'])) + '\n')
        print('Variance gravy: '+str(variance(proteins['Gravy'])))
    return True

def write_df(dataframe):
    dataframe.to_csv('../outputs/protein_db.csv', index_label='Protein Name')
    return True

if __name__ == '__main__':
    main()
