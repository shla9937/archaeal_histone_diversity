#!/usr/bin/env python3

import sys
import argparse
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import plotly

def main():
    parser = argparse.ArgumentParser(description='Import proteins df, cluster using DBSCAN, output df with clusters.')
    parser.add_argument('--df', required=True, type=str,
                        help='Pandas df in csv format that doesn\'t already contain cluster annotation.')
    parser.add_argument('--sil', required=True, type=str,
                        help='Silhouette score used for optimization of DBSCAN.')
    args = parser.parse_args()
    imports = pd.read_csv(args.df)
    imports.set_index('Protein Name', inplace=True)

    global features, cud_palette
    features = ['New Length', 'New pI', 'New Instability', 'New Gravy']
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"] 

    proteins = cluster(float(args.sil), imports)
    fig = plot(proteins)
    plotly.offline.plot(fig, filename='../outputs/dbscan_optimized_plot.html', auto_open=True)
    write_df(proteins, args.df)

def update_counter(final_sil, final_eps, final_min):
    sys.stdout.write("\rCurrent sil: {}, eps: {}, min: {} ".format(final_sil, final_eps, final_min))
    sys.stdout.flush()

def optimize_parameters(X, sil_target):
    final_eps = 0
    final_min = 0
    final_sil = 1
    for eps in np.arange(0.1, 2.0, 0.1):
        for min_samples in range(5, 100, 5):    
            db = DBSCAN(eps=eps, min_samples=int(min_samples)).fit(X)
            unique_clusters = np.unique(db.labels_)
            if len(unique_clusters) > 2 and len(unique_clusters) < 11:
                silhouette_avg = silhouette_score(X, db.labels_)
                if abs(silhouette_avg - sil_target) < abs(final_sil - sil_target):
                    update_counter(final_sil, final_eps, final_min)
                    final_eps = eps
                    final_min = min_samples
                    final_sil = silhouette_avg 
    update_counter(final_sil, final_eps, final_min)
    return final_eps, final_min

def cluster(sil_target, proteins):
    X = np.c_[proteins['New Length'], proteins['New pI'], proteins['New Instability'], proteins['New Gravy']]
    eps, min_samples = optimize_parameters(X, sil_target)
    db = DBSCAN(eps=eps, min_samples=int(min_samples)).fit(X)
    silhouette_avg = silhouette_score(X, db.labels_)
    print("\n Final sil: "+str(silhouette_avg)+", eps: "+str(eps)+", min: "+str(min_samples))
    print("Clusters: ", len(np.unique(db.labels_)) - 1)   
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    proteins['Cluster'] = db.labels_.astype(str)
    proteins.sort_values('Cluster', axis=0, ascending=True, inplace=True, kind='quicksort')
    return proteins

def plot(proteins):
    fig = px.scatter_3d(proteins, x='Length', y='pI', z='Gravy', color='Cluster', color_discrete_sequence=cud_palette, 
                        hover_name=proteins.index)
    fig.update_traces(marker=dict(size=3, opacity=0.6))
    fig.update_layout(title={'text': 'DBSCAN optimized for sil score'}, legend={'itemsizing': 'constant'}) 
    fig.update_scenes(camera_eye_x=-1.5)
    fig.update_scenes(camera_eye_y=-1.75)
    fig.update_scenes(camera_eye_z=0.15)
    fig.update_layout(width=1000, height=750)
    return fig

def write_df(dataframe, old_file):
    dataframe.to_csv(old_file, index_label='Protein Name')
    return True

if __name__ == '__main__':
    main()