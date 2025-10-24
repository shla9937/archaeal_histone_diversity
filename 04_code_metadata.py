#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import ast

def main():
    parser = argparse.ArgumentParser(description='Process protein data and plot interactions.')
    parser.add_argument('--df', type=str, required=True, help='Path to the CSV file containing protein df')
    parser.add_argument('--meta', type=str, required=True, help='Path to TSV containing metadata')
    args = parser.parse_args()
    
    global clust_features
    clust_features = ['New Length', 'New pI', 'New Instability', 'New Gravy']

    proteins = pd.read_csv(args.df)
    proteins.set_index('Protein Name', inplace=True)
    proteins['Cluster'] = proteins['Cluster'].astype(str)
    metadata = pd.read_csv(args.meta, sep='\t')
    locations = code_locations(metadata)
    output_locations(locations)
    proteins_by_species, nulls = merge_dfs(proteins, locations)
    proteins_by_species = add_strategy(proteins_by_species)
    output_species(proteins_by_species, nulls)
    write_strategies(proteins_by_species, proteins)

def code_locations(metadata):
    metadata.set_index('accession', inplace=True)
    rep_genomes = metadata['gtdb_genome_representative'].unique()
    locations = metadata.loc[metadata.index.isin(rep_genomes), 
                             ['coding_density', 'gc_percentage', 'genome_size', 'gtdb_taxonomy', 'ncbi_isolation_source']]
    locations['standard_location'] = locations['ncbi_isolation_source'].apply(determine_location)
    locations['pressure'] = locations['standard_location'].apply(determine_pressures)
    locations = add_tax(locations)
    locations = locations.drop(columns=['gtdb_taxonomy'])
    return locations

def add_tax(df):
    for idx, row in df.iterrows():
        if pd.notna(row['gtdb_taxonomy']):
            tax_levels = row['gtdb_taxonomy'].split(';')
            df.at[idx, 'Domain'] = tax_levels[0].lstrip('d__')
            df.at[idx, 'Phylum'] = tax_levels[1].lstrip('p__')
            df.at[idx, 'Class'] = tax_levels[2].lstrip('c__')
            df.at[idx, 'Order'] = tax_levels[3].lstrip('o__')
            df.at[idx, 'Family'] = tax_levels[4].lstrip('f__')
            df.at[idx, 'Genus'] = tax_levels[5].lstrip('g__')
            df.at[idx, 'Species'] = tax_levels[6].lstrip('s__')
    return df

def determine_location(location):
    if pd.isna(location):  # Handles nan values
        return ["none"]
    location = location.lower()
    locations = []
    
    if "acid" in location or "ph 4" in location or "ph 3" in location:
        locations.append("acid")
    if "soda" in location or "alkaline" in location:
        locations.append("alkaline")
    if "seep" in location:
        locations.append("cold seep")
    if "permafrost" in location:
        locations.append("permafrost")
    if "salt" in location or "brine" in location or "saline" in location:
        locations.append("salt")
    if "hot spring" in location or "hot pool" in location or "volcano" in location or "hot solfataric spring" in location or "geothermal" in location:
        locations.append("hot spring") 
    if "mine" in location or "mining" in location:
        locations.append("mine")  
    if "gut" in location or "intestine" in location or "rumen" in location or "ruminant" in location or "caecal" in location or "sheep" in location or "intestinal" in location or "caecum" in location:
        locations.append("gut")
    if "vent" in location or "hydrothermal" in location or "chimney" in location or "black smoker" in location:
        locations.append("hydrothermal vent")
    if "nuclear" in location or "olkiluoto" in location:
        locations.append("nuclear")
    if "sediment" in location or "mud" in location:
        locations.append("sediment")
    if "compost" in location:
        locations.append("compost")
    if "anaerobic digester" in location or "bioreactor" in location or "digestion" in location or "biodigester" in location or "reactor" in location or "fermentation" in location:
        locations.append("bioreactor/digester")
    if "waste" in location or "wwtp" in location or "sludge" in location:
        locations.append("waste")
    if "landfill" in location:
        locations.append("landfill")
    if "feces" in location or "fecal" in location or "stool" in location or "faeces" in location or "manure" in location:
        locations.append("feces")
    if "mouth" in location or "oral" in location or "dental" in location:
        locations.append("oral")
    if "estuary" in location:
        locations.append("estuary")
    if "ocean" in location or "marine" in location or "sea" in location or "atlantic" in location or "intertidal" in location or "tide-pools" in location or " bay" in location or "deep" in location or "depth" in location or "coastal" in location or "pacific" in location:
        locations.append("marine")
    if "soil" in location or "grass" in location:
        locations.append("soil")
    if "groundwater" in location or "well" in location or "aquifer" in location or "spring water" in location:
        locations.append("groundwater")
    if "marine" not in locations and ("water" in location or "aquatic" in location or "wetland" in location or "river" in location or "pond" in location or "peatland" in location or "lake" in location):
        locations.append("water")
    if "biofilm" in location:
        locations.append("biofilm")
    if ("oil" in location or "petro" in location) and "soil" not in location:
        locations.append("oil")
    if "lab" in location:
        locations.append("lab")
    if "frac" in location:
        locations.append("fracking")
    if "rock" in location or "aspo hrl" in location:
        locations.append("rock")
    elif len(locations) == 0 and ("none" in location or "metagenome" in location or location.strip() == ""):
        locations.append("none")
    elif len(locations) == 0:
        locations.append("other")
    return locations

def determine_pressures(locations):
    pressures = []

    for loc in locations:
        if loc == 'acid':
            pressures.extend(['acidic'])
        elif loc == 'alkaline':
            pressures.extend(['alkaline'])
        elif loc == 'cold seep':
            pressures.extend(['hyperbaric', 'cold', 'anaerobic'])
        elif loc == 'permafrost':
            pressures.extend(['cold'])
        elif loc == 'salt':
            pressures.extend(['saline'])
        elif loc == 'hot spring':
            pressures.extend(['hot'])
        elif loc == 'mine':
            pressures.extend(['contaminated', 'acidic'])
        elif loc == 'gut':
            pressures.extend(['acidic','competative'])
        elif loc == 'hydrothermal vent':
            pressures.extend(['hyperbaric', 'hot', 'anaerobic'])
        elif loc == 'nuclear':
            pressures.extend(['radioactive', 'contaminated'])
        elif loc == 'sediment':
            pressures.extend(['sediment'])
        elif loc == 'compost':
            pressures.extend(['nutrient-rich'])
        elif loc == 'bioreactor/digester':
            pressures.extend(['nutrient-rich', 'anaerobic'])
        elif loc == 'waste':
            pressures.extend(['nutrient-rich', 'anaerobic', 'contaminated'])
        elif loc == 'landfill':
            pressures.extend(['nutrient-rich', 'anaerobic', 'contaminated'])
        elif loc == 'feces':
            pressures.extend(['nutrient-rich', 'anaerobic'])
        elif loc == 'oral':
            pressures.extend(['competative'])
        elif loc == 'estuary':
            pressures.extend(['marine', 'freshwater'])
        elif loc == 'marine':
            pressures.extend(['marine'])
        elif loc == 'soil':
            pressures.extend(['soil'])
        elif loc == 'groundwater':
            pressures.extend(['freshwater'])
        elif loc == 'water':
            pressures.extend(['freshwater'])
        elif loc == 'biofilm':
            pressures.extend(['competative'])
        elif loc == 'oil':
            pressures.extend(['oil'])
        elif loc == 'lab':
            pressures.extend(['lab'])
        elif loc == 'fracking':
            pressures.extend(['oil'])
        elif loc == 'rock':
            pressures.extend(['nutrient-poor'])
        elif loc == 'none':
            pressures.extend(['none'])
        elif loc == 'other':
            pressures.extend(['other'])
    pressures = list(set(pressures))
    return pressures
    
def output_locations(locations):
    output_file_path = '../outputs/location_by_species.csv'
    locations.to_csv(output_file_path, index=True)
    return True   
     
def merge_dfs(proteins, locations):
    locations['Number of histones'] = 0
    locations['Clusters'] = [[] for _ in range(len(locations))]
    locations['Histones'] = [[] for _ in range(len(locations))]
    columns_order = ['Number of histones', 'Clusters', 'Histones'] + \
                    [col for col in locations.columns if col not in ['Number of histones', 'Clusters', 'Histones']]
    locations = locations[columns_order]
    nulls_by_species = pd.DataFrame(columns=['Number of proteins', 'Proteins', 'Species ID'])
    for index, row in proteins.iterrows():
        species = row['Species ID']
        cluster = row['Cluster']
        protein = index
        if cluster != '-1':
            if species in locations.index:
                locations.at[species, 'Number of histones'] += 1
                locations.at[species, 'Clusters'].append(cluster)
                locations.at[species, 'Histones'].append(protein)
            else:
                locations.loc[species] = {'Number of histones': 1,
                                          'Clusters': [cluster],
                                          'Histones': [protein]}
        else:
            if species in nulls_by_species.index:
                nulls_by_species.at[species, 'Number of proteins'] += 1
                nulls_by_species.at[species, 'Proteins'].append(protein)
            else:
                nulls_by_species.loc[species] = {'Number of proteins': 1,
                                                 'Proteins': [protein],
                                                 'Species ID': species}              
    return locations, nulls_by_species

def add_strategy(proteins_by_species):
    proteins_by_species['Clusters'] = proteins_by_species['Clusters'].apply(
        lambda x: x if isinstance(x, list) else [])
    proteins_by_species['Strategy'] = proteins_by_species['Clusters'].apply(generate_strategy)
    strategies = proteins_by_species['Strategy'].unique()
    singles = sorted([s for s in strategies if s.startswith('Single')])
    multiples = sorted([s for s in strategies if s.startswith('Multiple')])
    combinations = sorted([s for s in strategies if s.startswith('Combination')])
    no_histones = [s for s in strategies if s == 'No histones']
    category_order = singles + multiples + combinations + no_histones
    proteins_by_species['Strategy'] = pd.Categorical(
        proteins_by_species['Strategy'], categories=category_order, ordered=True)
    columns = list(proteins_by_species.columns)
    columns.insert(columns.index('Clusters'), columns.pop(columns.index('Strategy')))
    proteins_by_species = proteins_by_species[columns]
    return proteins_by_species

def generate_strategy(clusters):
    if not clusters:  
        return "No histones" 
    unique_clusters = sorted(set(clusters))
    if len(unique_clusters) == 1:
        if clusters.count(unique_clusters[0]) > 1:
            return f"Multiple {unique_clusters[0]}"
        else:
            return f"Single {unique_clusters[0]}"    
    elif len(unique_clusters) > 1:
        return "Combination " + ' + '.join(unique_clusters)

def output_species(proteins_by_species, nulls):
    nulls.to_csv('../outputs/null_species.csv', index=True)
    proteins_by_species.to_csv('../outputs/proteins_clustered_taxonomy_strategy_db.csv', index=True)
    return True 

def write_strategies(strategies_df, sequences_df):
    closest_proteins = {}
    strategy_dict = {}
 
    for _, row in strategies_df.iterrows():
        strategy = row['Strategy'].replace(" ", "_")
        clusters = row['Clusters']
        histones = row['Histones']
        for cluster in clusters:
            cluster_key = f"{strategy}_{cluster}"
            if cluster_key not in strategy_dict:
                strategy_dict[cluster_key] = []
    for _, row in strategies_df.iterrows():
        strategy = row['Strategy'].replace(" ", "_")
        clusters = row['Clusters']
        histones = row['Histones']
        for i, cluster in enumerate(clusters):
            cluster_key = f"{strategy}_{cluster}"
            strategy_dict[cluster_key].append(histones[i])
    for cluster_key, histone_list in strategy_dict.items():
        with open(f'../outputs/{cluster_key}.fa', 'w') as f:
            for protein_name in histone_list:
                if protein_name in sequences_df.index:
                    sequence = sequences_df.loc[protein_name, 'Sequence']
                    f.write(f'{protein_name}\n{sequence}\n')
    with open('../outputs/closest_to_com.fa', 'w') as f_com:
        for strategy, histone_list in strategy_dict.items():
            subset = sequences_df.loc[histone_list].dropna(subset=clust_features)
            if not subset.empty:
                com = subset[clust_features].mean().values
                distances = ((subset[clust_features] - com) ** 2).sum(axis=1)
                closest_idx = distances.idxmin()
                closest_sequence = sequences_df.loc[closest_idx, 'Sequence']
                closest_proteins[strategy] = closest_sequence
                f_com.write(f'{closest_idx}_{strategy}\n{closest_sequence}\n')
    return True

if __name__ == '__main__':
    main()

