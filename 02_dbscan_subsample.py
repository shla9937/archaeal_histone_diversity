#!/usr/bin/env python3

import sys
import argparse
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.offline as pyo
import statistics
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
    
    global clust_features, plot_features, cud_palette
    clust_features = ['New Length', 'New pI', 'New Instability', 'New Gravy']
    plot_features = ['Length', 'pI', 'Gravy']
    cud_palette = ["#999999","#0072B2","#56B4E9","#E69F00","#F0E442","#009E73","#D55E00","#CC79A7","#000000"]
    
    e, min_samples, bounds = get_bounds(imports, args.sil, args.n)
    final_bounds = test_bounds(imports, args.n, e, min_samples, bounds)
    if final_bounds is not None:
        proteins = assign_final_labels(imports, final_bounds)
        fig = plot(proteins)
        pyo.offline.plot(fig, filename=f'../outputs/dbscan_optimized_plot.html', auto_open=True)
        write_df(proteins, '../outputs/proteins_clustered_db.csv')
    else:
        print("Bounds were not found, try a different sample size")

def optimize_parameters(X, sil_target):
    final_eps = 0
    final_min = 0
    final_sil = 1
    for eps in np.arange(0.1, 2.0, 0.1):
        for min_samples in range(5, 100, 5):  
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
    db = DBSCAN(eps=e, min_samples=int(min_samples)).fit(sample_df[clust_features])
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
    all_bounds = {}
    initial = True
    print(f"Finding e and min")

    for run in range(3):
        sample_df = imports.sample(n=n_samples)
        X = sample_df[clust_features].values
        
        if initial:
            e, min_samples = optimize_parameters(X, sil_target)
            initial = False
        
        sample_df = get_clusters(sample_df, X, e, min_samples)
        fig = plot(sample_df)
        pyo.plot(fig, filename=f'../outputs/run_{run}_plot.html', auto_open=True)
        
        run_bounds = {}
        for cluster in sample_df['Cluster'].unique():
            if cluster == '-1':
                continue
            
            cluster_data = sample_df[sample_df['Cluster'] == cluster]
            cluster_bounds = {}
            for feature in clust_features:
                cluster_bounds[feature] = {
                    'min': cluster_data[feature].min(),
                    'max': cluster_data[feature].max(),
                    'center': (cluster_data[feature].min() + cluster_data[feature].max()) / 2}
            run_bounds[cluster] = cluster_bounds
        
        if not bounds:
            bounds = run_bounds
            all_bounds = {cluster: {feature: {'min': [], 'max': []} for feature in clust_features} for cluster in bounds}
        else:
            run_clusters = set(run_bounds.keys())
            existing_clusters = set(bounds.keys())
            
            if run_clusters != existing_clusters:
                sys.exit(f"Uneven number of clusters found: Run clusters = {run_clusters}, Existing clusters = {existing_clusters}.")

            for run_cluster, run_cluster_bounds in run_bounds.items():
                if run_cluster == '-1':
                    continue
                
                matching_cluster = None
                for existing_cluster, existing_bounds in bounds.items():
                    centroid = {feature: (run_cluster_bounds[feature]['min'] + run_cluster_bounds[feature]['max']) / 2.0 for feature in clust_features}
                    fits_all = True
                    for feature in clust_features:
                        if not (existing_bounds[feature]['min'] <= centroid[feature] <= existing_bounds[feature]['max']):
                            fits_all = False
                            break
                    if fits_all:
                        matching_cluster = existing_cluster
                        break
                
                if matching_cluster is not None:
                    for feature in clust_features:
                        feature_min = run_cluster_bounds[feature]['min']
                        feature_max = run_cluster_bounds[feature]['max']
                        all_bounds[matching_cluster][feature]['min'].append(feature_min)
                        all_bounds[matching_cluster][feature]['max'].append(feature_max)
                else:
                    print(f"Cluster {run_cluster} does not match any existing cluster.")

    final_bounds = {}
    for cluster, features_bounds in all_bounds.items():
        final_bounds[cluster] = {}
        for feature, bounds_values in features_bounds.items():
            if bounds_values['min']:
                final_bounds[cluster][feature] = {
                    'min': sum(bounds_values['min']) / len(bounds_values['min']),
                    'max': sum(bounds_values['max']) / len(bounds_values['max'])}
            else:
                final_bounds[cluster][feature] = {
                    'min': None,
                    'max': None}
    
    return e, min_samples, final_bounds
        
def test_bounds(imports, n_samples, e, min_samples, bounds):
    final_bounds = None
    passes = 0
    fails = 0
    for i in range(3):
        print(f"Running test {i + 1}/20")
        sample_df = imports.sample(n=n_samples)
        X = sample_df[clust_features].values
        sample_df = get_clusters(sample_df, X, e, min_samples)
        test_df = sample_df.copy()
        test_df['Cluster'] = sample_df.apply(lambda row: classify_point(row, bounds), axis=1)
        test_df = clean_df(test_df)
        plot_test = test_df.copy()
        test_df = test_df.reindex(sample_df.index)

        test_df['Cluster'] = test_df['Cluster'].astype(str)
        test_df['Cluster'] = pd.Categorical(test_df['Cluster'], categories=sample_df['Cluster'].cat.categories, ordered=True)

        same_cluster = (sample_df['Cluster'] == test_df['Cluster'])
        print(same_cluster.mean())
        if same_cluster.mean() > 0.95:
            passes += 1
        else:
            fails += 1

    print("Percent of passing tests: "+str((passes/(passes+fails))*100)+"%")
    
    if (passes/(passes+fails)) > 0.67:
        final_bounds = bounds

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
    fig = px.scatter_3d(proteins, x=plot_features[0], y=plot_features[1], z=plot_features[2], color='Cluster', 
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

def assign_final_labels(imports, final_bounds):
    imports['Cluster'] = imports.apply(lambda row: classify_point(row, final_bounds), axis=1)
    proteins = clean_df(imports)
    write_bounds_to_file(final_bounds, proteins)
    return proteins

def write_bounds_to_file(final_bounds, proteins):
    rows = []
    for cluster, bounds in final_bounds.items():
        row = {'Cluster': cluster}
        for feature, feature_bounds in bounds.items():
            original_feature_name = feature.replace('New ', '')
            original_data = proteins[original_feature_name].tolist()
            unstandardized_min = unstandardize_data(original_data, [feature_bounds['min']])[0]
            unstandardized_max = unstandardize_data(original_data, [feature_bounds['max']])[0]
            mean_value = proteins[proteins['Cluster'] == cluster][original_feature_name].mean()
            median_value = proteins[proteins['Cluster'] == cluster][original_feature_name].median()
            std_dev_value = proteins[proteins['Cluster'] == cluster][original_feature_name].std()
            row[f'{original_feature_name}_min'] = unstandardized_min
            row[f'{original_feature_name}_max'] = unstandardized_max
            row[f'{original_feature_name}_mean'] = mean_value
            row[f'{original_feature_name}_median'] = median_value
            row[f'{original_feature_name}_std_dev'] = std_dev_value
        rows.append(row)
    bounds_df = pd.DataFrame(rows)
    bounds_df.to_csv('../outputs/dbscan_bounds.csv', index=False)
    return True

def unstandardize_data(data_list, coords):
    unstandard_list = []
    list_stdev = statistics.stdev(data_list)
    list_mean = statistics.mean(data_list)
    for i in coords:
        new_i = (i * list_stdev) + list_mean
        unstandard_list.append(new_i)
    return unstandard_list

def write_df(dataframe, file_name):
    dataframe.to_csv(file_name, index_label='Protein Name')
    return True

if __name__ == '__main__':
    main()
