#!/usr/bin/env python3

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Aggregate taxa and strategies from input CSV")
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-c', '--cutoff', type=int, required=True, help='Cutoff for the minimum number of occurrences of each strategy')
    parser.add_argument('-t', '--taxonomy', required=True, help='Taxonomy file')
    
    args = parser.parse_args()
    
    aggregated_data = aggregate_taxa_strategies(args.input, args.cutoff)
    tax_dict = load_taxonomy(args.taxonomy)
    all_possible_nodes = create_all_possible_nodes(tax_dict)
    
    write_singles(aggregated_data, all_possible_nodes)
    
def aggregate_taxa_strategies(input_csv, cutoff):
    df = pd.read_csv(input_csv)
    taxonomic_columns = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'accession']
    strategy_column = 'Strategy'
    strategy_counts = df[strategy_column].value_counts()
    
    valid_strategies = strategy_counts[strategy_counts >= cutoff].index.tolist()
    print(f"Valid strategies (count >= {cutoff}): {valid_strategies}")
    
    aggregated_data = []
    
    for taxon in taxonomic_columns:
        unique_taxa = df[taxon].unique()
        for value in unique_taxa:
            if taxon == 'accession':
                row = {'accession': value}
            else:
                row = {taxon: f"{taxon[0].lower()}__{value}"}
            
            for strategy in valid_strategies:
                associated = ((df[taxon] == value) & (df[strategy_column] == strategy)).any()
                row[strategy] = 1 if associated else 0 
            aggregated_data.append(row)
    return aggregated_data

def load_taxonomy(tax_file):
    tax_dict = {}
    with open(tax_file, 'r') as tax_file:
        for line in tax_file:
            line_list = line.strip().split('\t')
            tax_dict[line_list[0]] = line_list[1].split(';')
    return tax_dict

def create_all_possible_nodes(tax_dict):
    all_nodes = set()

    for specium in tax_dict:
        clades = tax_dict[specium]
        for i in range(len(clades)):
            # Iterate over all possible combinations of clades
            for j in range(i + 1, len(clades) + 1):
                # Create node name for the clade
                node_name = '|'.join(clades[i:j])
                all_nodes.add(node_name)

    return sorted(all_nodes)  # Return a sorted list of all unique nodes

def write_singles(aggregated_data, all_possible_nodes):
    # Define the header for the iTOL file
    itol_header = """DATASET_BINARY
SEPARATOR COMMA
DATASET_LABEL,Singles
COLOR,#aaaaaa
FIELD_COLORS,#56B4E9,#D55E00,#E69F00,#009E73
FIELD_LABELS,1,2,3,5
FIELD_SHAPES,2,2,2,2
DASHED_LINES,1
DATA"""

    with open('itol_singles.txt', 'w') as f:
        f.write(itol_header + '\n')
        
        node_dict = {node: [0] * len(aggregated_data[0].keys()) for node in all_possible_nodes}
        
        for entry in aggregated_data:
            taxon_label = list(entry.values())[0]
            strategy_values = [entry[strategy] for strategy in sorted(entry.keys()) if 'Single' in strategy]
            for node in all_possible_nodes:
                # Check if the taxon_label is part of the node
                if taxon_label in node:  # Adjusted from taxon_label to use `in` for more flexibility
                    node_index = list(entry.values())[0]  # Use the first key or appropriate method to get the correct index
                    node_dict[node] = [x + y for x, y in zip(node_dict[node], strategy_values)]

        # Write nodes and their strategy values to the file
        for node, values in node_dict.items():
            values_str = ', '.join(map(str, values))
            f.write(f"{node}, {values_str}\n")


if __name__ == "__main__":
    main()
