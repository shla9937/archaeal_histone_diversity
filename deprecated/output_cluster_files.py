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
    cluster_data = cluster_proteins(imported_proteins_atts)
    write_clusters(cluster_data)
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

def cluster_proteins(imported_proteins_atts):
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
    
    cluster_data = []
    cluster_null = []
    cluster_0 = []
    cluster_1 = []
    cluster_2 = []
    cluster_3 = []
    cluster_4 = []
    cluster_5 = []
    cluster_6 = []
    
    for label, protein in zip(labels, imported_proteins_atts):
        new_protein = []
        new_protein.append(protein[0])
        new_protein.append(protein[1])
        if label == -1:
            cluster_null.append(new_protein)
        elif label == 0:
            cluster_0.append(new_protein)
        elif label == 1:
            cluster_1.append(new_protein)
        elif label == 2:
            cluster_2.append(new_protein)
        elif label == 3:
            cluster_3.append(new_protein)
        elif label == 4:
            cluster_4.append(new_protein)
        elif label == 5:
            cluster_5.append(new_protein)
        elif label == 6:
            cluster_6.append(new_protein)
            
    cluster_data.append(cluster_0)
    cluster_data.append(cluster_1)
    cluster_data.append(cluster_2)
    cluster_data.append(cluster_3)
    cluster_data.append(cluster_4)
    cluster_data.append(cluster_5)
    cluster_data.append(cluster_6)
    cluster_data.append(cluster_null)
        
    return(cluster_data)

def write_clusters(cluster_data):
    cluster_num = 0
    for cluster in cluster_data:
        if cluster_num == len(cluster_data)-1:
            f = open('cluster_null.fa', 'a')
        else:
            f = open('cluster_{}.fa'.format(cluster_num), 'a')
        for protein in cluster:
            f.write(protein[0]+'\n'+protein[1]+'\n')
        f.close()
        cluster_num += 1
    return True


if __name__ == '__main__':
    main()