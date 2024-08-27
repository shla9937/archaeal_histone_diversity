import sys
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
import matplotlib 
from matplotlib import pyplot as plt
import statistics
from statistics import stdev
from statistics import mean


def main():
    input_file = sys.argv[1]
    standard = sys.argv[2]
    imported_proteins = import_file(input_file)
    imported_proteins_atts = calc_atts(imported_proteins)
    if standard == 'False':
        plot_proteins(imported_proteins_atts)
    elif standard == 'True':
        plot_standard_proteins(imported_proteins_atts)
    return True

def import_file(input_file):
    imported_proteins = []
    sequence = None
    seq_len = 0
    gene_count = 0

    with open(input_file, 'r') as fastafile:
        for line in fastafile:
            if '>' in line:
                gene_count += 1
                if sequence is not None:
                    seq_len = len(sequence)
                    protein.append(sequence)
                    protein.append(seq_len)
                    sequence = None
                    imported_proteins.append(protein)
                protein = []
                protein.append(line.rstrip('\n'),)
            elif (len(line) != 0) and (sequence is not None):
                sequence += (line.rstrip('\n'))
            elif (len(line) != 0) and (sequence is None):
                sequence = (line.rstrip('\n'))
        seq_len = len(sequence)
        protein.append(sequence)
        protein.append(seq_len)
        imported_proteins.append(protein)
    return imported_proteins

def calc_atts(imported_proteins):
    imported_proteins_atts = []
    for protein in imported_proteins:
        new_protein = []
        X = ProtParam.ProteinAnalysis(protein[1])
        new_protein.append(protein[0])
        new_protein.append(protein[1])
        new_protein.append(protein[2]) 
        new_protein.append(X.aromaticity())
        new_protein.append(X.isoelectric_point()) 
        new_protein.append(X.secondary_structure_fraction()[0])
        imported_proteins_atts.append(new_protein)
    return imported_proteins_atts

def plot_proteins(imported_proteins_atts):
    fig, ax = plt.subplots(1, 3)
    isos = []
    lengths = []
    helices = []
    aromatics = []
    for protein in imported_proteins_atts:
        isos.append(protein[4])
        lengths.append(protein[2])
        aromatics.append(protein[3])
        helices.append(protein[-1])    
    ax[0].scatter(lengths, isos)
    ax[0].set_title('Len vs. pI')
    ax[1].scatter(lengths, helices)
    ax[1].set_title('Len vs. helices')
    ax[2].scatter(lengths, aromatics)
    ax[2].set_title('Len vs. aromatics')
    plt.tight_layout()
    plt.show()

def standardize_data(data_list):
    standard_list = []
    list_stdev = statistics.stdev(data_list)
    list_mean = statistics.mean(data_list)
    for i in data_list:
        new_i = (i - list_mean) / list_stdev
        standard_list.append(new_i)
    return standard_list

def plot_standard_proteins(imported_proteins_atts):
    fig2, ax = plt.subplots(1, 3)
    isos = []
    lengths = []
    helices = []
    aromatics = []
    for protein in imported_proteins_atts:
        isos.append(protein[4])
        lengths.append(protein[2])
        aromatics.append(protein[3])
        helices.append(protein[-1])
    lengths = standardize_data(lengths)
    helices = standardize_data(helices)
    aromatics = standardize_data(aromatics)
    isos = standardize_data(isos)
    ax[0].scatter(lengths, isos)
    ax[0].set_title('Len vs. pI')
    ax[1].scatter(lengths, helices)
    ax[1].set_title('Len vs. helices')
    ax[2].scatter(lengths, aromatics)
    ax[2].set_title('Len vs. aromatics')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
