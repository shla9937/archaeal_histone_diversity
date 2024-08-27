#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import statistics
import plotly
import plotly.express as px
import plotly.graph_objects as go

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process protein data and plot interactions.')
    parser.add_argument('--df', type=str, required=True, help='Path to the CSV file containing protein df')
    parser.add_argument('--tax', type=str, required=True, help='Path to the taxonomy file')
    args = parser.parse_args()
    
    # Load data
    proteins = pd.read_csv(args.df)
    proteins.set_index('Protein Name', inplace=True)
    proteins['Cluster'] = proteins['Cluster'].astype(str)
    # Process data
    proteins = add_tax(proteins, args.tax)
    interaction_dict, pairs_dict = find_interactions(proteins)
    
    # Plot network
    color_palette()
    fig = plot(proteins)
    fig = plot_network(fig, interaction_dict, pairs_dict, proteins, True)
    fig, cluster_centers = add_centers(proteins, fig)
    fig, centroids = add_centroids(proteins, cluster_centers, fig)
    fig = add_dropdown(proteins, fig)
    fig.show()
    plotly.offline.plot(fig, filename='../outputs/3d_interactive_plot.html', auto_open=True)
    write_clusters(proteins)
    write_df(proteins, args.df)
    
def add_tax(proteins, tax_file):
    tax_dict = {}
    with open(tax_file, 'r') as tax_file:
        for line in tax_file:
            line_list = line.strip().split('\t')
            tax_dict[line_list[0]] = line_list[1].split(';')
    for protein in proteins.index:
        species = proteins.loc[protein, 'Species ID']
        for idx, level in enumerate(['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']):
            proteins.at[protein, level] = tax_dict[species][idx].lstrip(f"{level[0].lower()}__")
    return proteins

def color_palette():
    global tab_20, tab_20_pastel
    old_tab_20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 179, 71),  
                  (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
                  (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
                  (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
                  (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    tab_20_pastel, tab_20 = [], []
    for idx, color in enumerate(old_tab_20):
        if idx % 2 == 1:
            tab_20_pastel.append('rgb'+str(color))
        else:
            tab_20.append('rgb'+str(color))
    return True

def unstandardize_data(data_list, coords):
    unstandard_list = []
    list_stdev = statistics.stdev(data_list)
    list_mean = statistics.mean(data_list)
    for i in coords:
        new_i = (i * list_stdev) + list_mean
        unstandard_list.append(new_i)
    return unstandard_list

def find_interactions(proteins):
    species = proteins['Species ID'].unique()
    interaction_dict = {}
    pairs_dict = {}

    protein_to_cluster = proteins['Cluster'].to_dict()

    for specium in species:
        interaction_dict[specium] = proteins[proteins['Species ID'] == specium].index.astype(str).tolist()
    
    for specium in interaction_dict:
        prime_list = interaction_dict[specium]
        query_list = prime_list.copy()  # Copy the list to avoid modification during iteration
        pairs = []
        
        for prime_protein in prime_list:
            for query_protein in query_list:
                if prime_protein != query_protein:
                    pairs.append([prime_protein, query_protein])
            query_list.remove(prime_protein)
        
        for pair in pairs:
            # Use the protein IDs to look up cluster values from the mapping
            clust_list = [protein_to_cluster.get(pair[0], None), protein_to_cluster.get(pair[1], None)]
            clust_set = set(clust_list)
            # Remove None values from clustering set
            clust_set.discard(None)
            clust_sort = sorted(clust_set)
            clust_tup = tuple(clust_sort)
            
            if clust_tup:
                if clust_tup in pairs_dict:
                    pairs_dict[clust_tup] += [pair]
                else: 
                    pairs_dict[clust_tup] = [pair]

    return interaction_dict, pairs_dict

def plot_network(fig, interaction_dict, pairs_dict, proteins, remove_null):
    inter_xs = []
    inter_ys = []
    inter_zs = []
    for interaction in sorted(pairs_dict):
        if remove_null is True and '-1' in interaction:
            continue 
        line_xs = []
        line_ys = []
        line_zs = []
        for pair in pairs_dict[interaction]:
            line_xs.append(proteins.loc[pair[0], 'Length']) 
            line_xs.append(proteins.loc[pair[1], 'Length'])
            line_xs.append(None)
            line_ys.append(proteins.loc[pair[0],'pI'])
            line_ys.append(proteins.loc[pair[1],'pI'])
            line_ys.append(None)
            line_zs.append(proteins.loc[pair[0], 'Gravy'])
            line_zs.append(proteins.loc[pair[1], 'Gravy'])
            line_zs.append(None)
        if len(interaction) == 1:
            cluster = str(interaction[0])
            trace_name = 'Cluster '+cluster+' interactions'
            color = tab_20_pastel[int(cluster)+1]
            fig.add_scatter3d(x=line_xs, y=line_ys, z=line_zs, mode='lines', 
                  line=dict(color=color, width=0.2), name=trace_name)
        else:
            for x, y, z in zip(line_xs, line_ys, line_zs):
                inter_xs.append(x)
                inter_ys.append(y)
                inter_zs.append(z)
    fig.add_scatter3d(x=inter_xs, y=inter_ys, z=inter_zs, mode='lines', 
                      line=dict(color='#aaaaaa', width=0.2), name='Inter-cluster interactions')
    
    return fig

def add_dropdown(proteins, fig):
    dropdown_list = []
    traces = list(fig.data)
    new_traces = proteins['Phylum'].unique().tolist()
    dropdown_list.append(dict(args=[{"visible": [True] + [False]*len(traces+new_traces)}], 
                         label="None", method="update"))
    idx = 0
    for phylum in proteins['Phylum'].unique():
        phyla_indexes = proteins[proteins['Phylum'] == phylum].index     
        trace = go.Scatter3d(
            x=proteins.loc[phyla_indexes, 'Length'],
            y=proteins.loc[phyla_indexes, 'pI'],
            z=proteins.loc[phyla_indexes, 'Gravy'],
            hovertext=phyla_indexes,
            mode='markers', hoverinfo=('text+x+y+z'),
            marker=dict(size=3, color='#000000'),
            name=phylum,
            visible=False)
        fig.add_traces(trace)
        dropdown_list.append(dict(args=[{"visible": [True]*len(traces) + idx*[False] + [True] + (len(traces+new_traces)-idx-1)*[False]}], 
                                  label=phylum, method="update"))
        idx += 1
    fig.update_layout(updatemenus=[dict(buttons=dropdown_list, x=0, y=0.95)])
    fig.update_layout(annotations=[dict(text="Phylum:", showarrow=False, x=-0.29, y=0.99)])
    return fig
    
def euclidean_distance(point1, point2):
    return np.sqrt(np.sum((point1 - point2)**2))
    
def add_centroids(proteins, cluster_centers, fig):
    columns_to_mean = ['New Length', 'New pI', 'New Gravy']
    centroids = {cluster: [] for cluster in cluster_centers}
    for cluster, center in cluster_centers.items():
        center = center[columns_to_mean]
        min_distance = float('inf')
        closest_centroid = None
        for protein, row in proteins.iterrows():
            point = row[columns_to_mean]
            distance = euclidean_distance(point, center)
            if distance < min_distance:
                min_distance = distance
                closest_centroid = protein
        centroids[cluster] = closest_centroid

    closest_x_values = []
    closest_y_values = []
    closest_z_values = []
    closest_proteins = []
    clusters = []
    f = open('../outputs/centroids.fa', 'w')
    for cluster, protein in centroids.items():
        closest_x_values.append(proteins.loc[protein, 'Length'])
        closest_y_values.append(proteins.loc[protein, 'pI'])
        closest_z_values.append(proteins.loc[protein, 'Gravy'])
        closest_proteins.append(protein)
        clusters.append(int(cluster))
        f.write(protein+'_centroid_of_cluster_'+str(cluster)+'\n'+proteins.loc[protein, 'Sequence']+'\n')
    f.close() 
    
    fig.add_scatter3d(x=closest_x_values, y=closest_y_values, z=closest_z_values, hovertext=closest_proteins, 
                      hoverinfo=('text+x+y+z'), marker=dict(color=tab_20_pastel[1:], size=10, opacity=0.6), 
                      mode='markers', name='Closest points')
    return fig, centroids

def add_centers(proteins, fig):
    cluster_points = {}
    columns_to_mean = ['New Length', 'New Aromaticity', 'New pI', 'New Helix', 'New Turn', 'New Sheet', 'New Instability', 'New Flexibility', 'New Gravy']
    for cluster in proteins['Cluster'].unique():
        if cluster == '-1':
            continue
        if cluster not in cluster_points:
            cluster_points[cluster] = []
        cluster_data = proteins[proteins['Cluster'] == cluster][columns_to_mean]
        cluster_points[cluster].append(cluster_data)  
    cluster_centers = {}
    for cluster in cluster_points:
        if cluster not in cluster_centers:
            combined_data = pd.concat(cluster_points[cluster], axis=0)
            cluster_centers[cluster] = combined_data.mean()
    x_values = unstandardize_data(proteins['Length'], [cluster_centers[cluster]['New Length'] for cluster in cluster_centers])
    y_values = unstandardize_data(proteins['pI'], [cluster_centers[cluster]['New pI'] for cluster in cluster_centers])
    z_values = unstandardize_data(proteins['Gravy'], [cluster_centers[cluster]['New Gravy'] for cluster in cluster_centers])
    
    fig.add_scatter3d(x=x_values, y=y_values, z=z_values, mode='markers', name='Centroids', marker=dict(color='#444444', size=6))
    return fig, cluster_centers

def plot(proteins):
    fig = px.scatter_3d(proteins, x='Length', y='pI', z='Gravy', color='Cluster', color_discrete_sequence=tab_20_pastel, 
                        hover_name=proteins.index, hover_data=['Species','Genus','Family','Order','Class','Phylum'])
    fig.update_traces(marker=dict(size=3, opacity=0.6))
    fig.update_layout(title={'text': 'Protein clustering'}, legend= {'itemsizing': 'constant'})
    fig.update_scenes(camera_eye_x=-1.5)
    fig.update_scenes(camera_eye_y=-1.75)
    fig.update_scenes(camera_eye_z=0.15)
    fig.update_layout(width=1000, height=750)
    return fig

def write_clusters(proteins):
    for cluster in proteins['Cluster'].unique():
        f = open('../outputs/cluster_{}.fa'.format(cluster), 'w')
        subset = proteins[proteins['Cluster'] == cluster]
        indexes = subset.index
        sequences = subset['Sequence']
        for index, sequence in zip(indexes, sequences):
            f.write(index+'\n'+sequence+'\n')
        f.close()
    
    g = open('../outputs/cluster_all_non_null.fa', 'w')
    subset = proteins[proteins['Cluster'] != str(-1)]
    indexes = subset.index
    sequences = subset['Sequence']
    for index, sequence in zip(indexes, sequences):
        g.write(index+'\n'+sequence+'\n')
    g.close()

    h = open('../outputs/cluster_all.fa', 'w')
    indexes = proteins.index
    sequences = proteins['Sequence']
    for index, sequence in zip(indexes, sequences):
        h.write(index+'\n'+sequence+'\n')
    h.close()    
    return True

def write_df(dataframe, old_file):
    dataframe.to_csv(old_file, index_label='Protein Name')
    return True

if __name__ == '__main__':
    main()

