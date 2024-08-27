#!/usr/bin/env python3

import sys
import argparse
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.offline as pyo
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def main():
    parser = argparse.ArgumentParser(description='Import proteins df, cluster using DBSCAN, output df with clusters.')
    parser.add_argument('--df', required=True, type=str,
                        help='Pandas df in csv format that doesn\'t already contain cluster annotation.')
    parser.add_argument('--sil', required=True, type=float,
                        help='Silhouette score used for optimization of DBSCAN.')
    parser.add_argument('--n', required=True, type=int,
                        help='Number of points to sample for each trial.')
    
    args = parser.parse_args()
    imports = pd.read_csv(args.df)
    imports.set_index('Protein Name', inplace=True)
    
    global features, cud_palette
    features = ['New Length', 'New pI', 'New Instability', 'New Gravy']
    cud_palette = ["#999999","#0072B2","#56B4E9","#E69F00","#F0E442","#009E73","#D55E00","#CC79A7","#000000"]
    
    e, min_samples, bounds = get_bounds(imports, args.sil, args.n)
    final_bounds = test_bounds(imports, args.n, e, min_samples, bounds)
    # if final_bounds not None:
    #     proteins = apply_bounds(imports, final_bounds)
    #     fig = plot(proteins)
    #     plotly.offline.plot(fig, filename='../outputs/dbscan_optimized_plot.html', auto_open=True)
    #     write_df(proteins, args.df)
    # else:
    #     print("Bounds were not found, try a different sample size")

def optimize_parameters(X, sil_target):
    final_eps = 0
    final_min = 0
    final_sil = 1
    # for eps in np.arange(0.1, 2.0, 0.1):
    #     for min_samples in range(5, 100, 5):  
    for eps in np.arange(0.2, 0.6, 0.1):
        for min_samples in range(20, 50, 10):  
            db = DBSCAN(eps=eps, min_samples=int(min_samples)).fit(X)
            unique_clusters = np.unique(db.labels_)
            num_clusters = len(unique_clusters)
            if len(unique_clusters) > 2:
                silhouette_avg = silhouette_score(X, db.labels_)
                if abs(silhouette_avg - sil_target) < abs(final_sil - sil_target):
                    update_counter(final_sil, final_eps, final_min, num_clusters)
                    final_eps = eps
                    final_min = min_samples
                    final_sil = silhouette_avg 
    update_counter(final_sil, final_eps, final_min, num_clusters)
    print(" ")
    return final_eps, final_min

def get_clusters(sample_df, X, e, min_samples):
    db = DBSCAN(eps=e, min_samples=int(min_samples)).fit(sample_df[features])
    sample_df['Cluster'] = db.labels_.astype(str)
    non_noise_df = sample_df[sample_df['Cluster'] != '-1']
    cluster_counts = non_noise_df['Cluster'].value_counts()
    sorted_clusters = cluster_counts.index
    cluster_mapping = {old_label: str(i + 1) for i, old_label in enumerate(sorted_clusters)}
    sample_df['Cluster'] = sample_df['Cluster'].map(cluster_mapping).fillna('-1')
    sample_df['Cluster'] = pd.Categorical(sample_df['Cluster'], categories=['-1'] + [str(i + 1) for i in range(len(sorted_clusters))], ordered=True)
    sample_df.sort_values('Cluster', axis=0, ascending=True, inplace=True, kind='quicksort')
    return sample_df

def get_bounds(imports, sil_target, n_samples):
    bounds = {}
    all_bounds = {}  # To store bounds from all runs for averaging
    initial = True
    print(f"Finding e and min")

    for run in range(3):
        sample_df = imports.sample(n=n_samples)
        X = sample_df[features].values
        
        if initial:
            e, min_samples = optimize_parameters(X, sil_target)
            initial = False
        
        sample_df = get_clusters(sample_df, X, e, min_samples)
        fig = plot(sample_df)
        pyo.plot(fig, filename=f'../outputs/run_{run}_plot.html', auto_open=True)
        
        run_bounds = {}
        for cluster in sample_df['Cluster'].unique():
            if cluster == '-1':
                continue  # Skip noise points
            
            cluster_data = sample_df[sample_df['Cluster'] == cluster]
            cluster_bounds = {}
            for feature in features:
                cluster_bounds[feature] = {
                    'min': cluster_data[feature].min(),
                    'max': cluster_data[feature].max(),
                    'center': (cluster_data[feature].min() + cluster_data[feature].max()) / 2}
            run_bounds[cluster] = cluster_bounds
        
        if not bounds:
            bounds = run_bounds
            # Initialize all_bounds with the first run's bounds
            all_bounds = {cluster: {feature: {'min': [], 'max': []} for feature in features} for cluster in bounds}
        else:
            # Compare clusters from this run to existing bounds
            run_clusters = set(run_bounds.keys())
            existing_clusters = set(bounds.keys())
            
            if run_clusters != existing_clusters:
                sys.exit(f"Uneven number of clusters found: Run clusters = {run_clusters}, Existing clusters = {existing_clusters}.")
            
            # Update existing bounds with new information
            for run_cluster, run_cluster_bounds in run_bounds.items():
                if run_cluster == '-1':
                    continue  # Skip noise points
                
                matching_cluster = None
                for existing_cluster, existing_bounds in bounds.items():
                    # Calculate centroid of the run_cluster
                    centroid = {feature: (run_cluster_bounds[feature]['min'] + run_cluster_bounds[feature]['max']) / 2.0 for feature in features}
                    
                    # Check if the centroid fits within the bounds of an existing cluster
                    fits_all = True
                    for feature in features:
                        if not (existing_bounds[feature]['min'] <= centroid[feature] <= existing_bounds[feature]['max']):
                            fits_all = False
                            break
                    if fits_all:
                        matching_cluster = existing_cluster
                        break
                
                if matching_cluster is not None:
                    # Store bounds for averaging
                    for feature in features:
                        feature_min = run_cluster_bounds[feature]['min']
                        feature_max = run_cluster_bounds[feature]['max']
                        all_bounds[matching_cluster][feature]['min'].append(feature_min)
                        all_bounds[matching_cluster][feature]['max'].append(feature_max)
                else:
                    print(f"Cluster {run_cluster} does not match any existing cluster.")
    
    # Calculate final bounds by averaging the stored min and max values
    final_bounds = {}
    for cluster, features_bounds in all_bounds.items():
        final_bounds[cluster] = {}
        for feature, bounds_values in features_bounds.items():
            if bounds_values['min']:  # Ensure there are values to average
                final_bounds[cluster][feature] = {
                    'min': sum(bounds_values['min']) / len(bounds_values['min']),
                    'max': sum(bounds_values['max']) / len(bounds_values['max'])}
            else:
                # If no values were collected, handle appropriately
                final_bounds[cluster][feature] = {
                    'min': None,
                    'max': None}
    
    return e, min_samples, final_bounds
        
def test_bounds(imports, n_samples, e, min_samples, bounds):
    final_bounds = None
    passes = 0
    fails = 0
    for i in range(20):
        print(f"Running test {i + 1}/20")
        sample_df = imports.sample(n=n_samples)
        X = sample_df[features].values
        sample_df = get_clusters(sample_df, X, e, min_samples)
        test_df = sample_df.copy()
        test_df['Cluster'] = sample_df.apply(lambda row: classify_point(row, bounds), axis=1)
        test_df = clean_df(test_df)
        plot_test = test_df.copy()
        test_df = test_df.reindex(sample_df.index)
        same_cluster = (sample_df['Cluster'] == test_df['Cluster'])
        print(same_cluster.mean())
        if same_cluster.mean() > 0.95:
            passes += 1
        else:
            fails += 1

    print("Percent of passing tests: "+str((passes/(passes+fails))*100)+"%")
    fig = plot(sample_df)
    pyo.plot(fig, filename=f'../outputs/sample_plot.html', auto_open=True)
    fig = plot(plot_test)
    pyo.plot(fig, filename=f'../outputs/test_plot.html', auto_open=True)
    
    return final_bounds

def clean_df(df):
    non_noise_df = df[df['Cluster'] != '-1']
    cluster_counts = non_noise_df['Cluster'].value_counts()
    sorted_clusters = cluster_counts.index
    cluster_mapping = {old_label: str(i + 1) for i, old_label in enumerate(sorted_clusters)}
    df['Cluster'] = df['Cluster'].map(cluster_mapping).fillna('-1')
    df['Cluster'] = pd.Categorical(df['Cluster'], categories=['-1'] + [str(i + 1) for i in range(len(sorted_clusters))], ordered=True)
    df.sort_values('Cluster', axis=0, ascending=True, inplace=True, kind='quicksort')
    return df

def classify_point(row, bounds):
    for cluster, cluster_bounds in bounds.items():
        fits_cluster = True
        for feature, feature_bounds in cluster_bounds.items():
            if not (feature_bounds['min'] <= row[feature] <= feature_bounds['max']):
                fits_cluster = False
                break
        if fits_cluster:
            return cluster
    return '-1'

def plot(proteins):
    fig = px.scatter_3d(proteins, x='Length', y='pI', z='Gravy', color='Cluster', 
                        color_discrete_sequence=cud_palette, category_orders={'Cluster': proteins['Cluster'].unique()}, 
                        hover_name=proteins.index)
    fig.update_traces(marker=dict(size=3, opacity=0.6))
    fig.update_layout(title={'text': 'DBSCAN plot'}, legend={'itemsizing': 'constant'}) 
    fig.update_scenes(camera_eye_x=-1.5)
    fig.update_scenes(camera_eye_y=-1.75)
    fig.update_scenes(camera_eye_z=0.15)
    fig.update_layout(width=1000, height=750)
    return fig

def update_counter(final_sil, final_eps, final_min, unique_clusters):
    sys.stdout.write(f"\rCurrent sil: {final_sil}, eps: {final_eps}, min: {final_min}, clusters: {unique_clusters}") 
    sys.stdout.flush()

def write_df(dataframe, old_file):
    dataframe.to_csv(old_file, index_label='Protein Name')
    return True

if __name__ == '__main__':
    main()












# def assign_final_labels(proteins, cluster_ranges, cluster_centers):
#     final_labels = []

#     for _, row in proteins.iterrows():
#         point = row[['New Length', 'New pI', 'New Instability', 'New Gravy']].values
#         assigned = False

#         for trial_num, ranges in cluster_ranges.items():
#             for cluster, cluster_range in ranges.items():
#                 if (cluster_range['New Length']['min'] <= point[0] <= cluster_range['New Length']['max'] and
#                     cluster_range['New pI']['min'] <= point[1] <= cluster_range['New pI']['max'] and
#                     cluster_range['New Instability']['min'] <= point[2] <= cluster_range['New Instability']['max'] and
#                     cluster_range['New Gravy']['min'] <= point[3] <= cluster_range['New Gravy']['max']):
#                     final_labels.append(str(cluster))  # Convert cluster number to string
#                     assigned = True
#                     break
#             if assigned:
#                 break

#         if not assigned:
#             final_labels.append('-1')  # Unassigned points

#     proteins['Final Cluster'] = final_labels

#     # Count occurrences of each cluster
#     cluster_counts = proteins['Final Cluster'].value_counts()
#     cluster_mapping = {cluster: str(rank + 1) for rank, (cluster, _) in enumerate(cluster_counts.items())}

#     # Map final cluster labels to most popular
#     proteins['Final Cluster'] = proteins['Final Cluster'].map(cluster_mapping).fillna('-1').astype(str)

#     # Plot final results using original feature values
#     final_fig = px.scatter_3d(proteins, x='Length', y='pI', z='Gravy', color='Final Cluster',
#                              color_discrete_sequence=px.colors.qualitative.Plotly, title='Final DBSCAN Clustering')
#     pyo.plot(final_fig, filename='../outputs/final_dbscan_labeled.html', auto_open=True)

#     # Save final labels to CSV
#     proteins.to_csv('../outputs/final_labeled_proteins.csv', index_label='Protein Name')
#     return proteins