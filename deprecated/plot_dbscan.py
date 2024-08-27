import sys
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
import matplotlib 
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import statistics
from statistics import stdev
from statistics import mean
from sklearn.cluster import DBSCAN
from sklearn import metrics
import numpy as np


def main():
    input_file = 'archaeal_histones_20191211.fa'
    imported_proteins = import_file(input_file)
    imported_proteins_atts = calc_atts(imported_proteins)
    plot_proteins(imported_proteins_atts)
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

def standardize_data(data_list):
    standard_list = []
    list_stdev = statistics.stdev(data_list)
    list_mean = statistics.mean(data_list)
    for i in data_list:
        new_i = (i - list_mean) / list_stdev
        standard_list.append(new_i)
    return standard_list

def plot_proteins(imported_proteins_atts):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
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
    
    X = np.c_[lengths, isos, helices]
    db = DBSCAN(eps=0.3, min_samples=10).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
        
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=labels, s=3, cmap='hsv')
    ax.set_xlabel('Lengths')
    ax.set_ylabel('pIs')
    ax.set_zlabel('Helices')
    ax.set_title('Archaeal histones')
    
    #plt.savefig('DBSCAN_3D.png', dpi=300)
    plt.show()
    
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    print("Silhouette Coefficient: %0.3f"
          % metrics.silhouette_score(X, labels))


if __name__ == '__main__':
    main()
