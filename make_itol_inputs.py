#!/usr/bin/env python3

import pandas as pd
import argparse
from ete3 import Tree

def main():
    parser = argparse.ArgumentParser(description="Aggregate taxa and strategies from input CSV")
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-c', '--cutoff', type=int, required=True, help='Cutoff for the minimum number of occurrences of each strategy')
    parser.add_argument('-t', '--tree', required=True, help='Tree file')
    args = parser.parse_args()
    
    aggregated_data = aggregate_taxa_strategies(args.input, args.cutoff)
    extra_node_ids = extract_node_ids(args.tree)
    write_singles(aggregated_data, extra_node_ids)
    write_multis(aggregated_data, extra_node_ids)
    write_combo12(aggregated_data, extra_node_ids)
    write_combo34(aggregated_data, extra_node_ids)
    
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
                row[strategy] = 1 if associated else -1 
            aggregated_data.append(row)
    return aggregated_data

def extract_node_ids(tree_file):
    tree = Tree(tree_file, format=1, quoted_node_names=True)
    extra_node_ids = []
    
    for node in tree.traverse():
        if node.name:
            if isinstance(node.name, str) and not node.name.isdigit() and ';' in node.name:
                extra_node_ids.append(node.name.replace('; ', '|'))

    cleaned_entries = []
    for entry in extra_node_ids:
        if ':' in entry:
            cleaned_entries.append(entry.split(':', 1)[1])
        else:
            cleaned_entries.append(entry)
    return cleaned_entries

def write_singles(aggregated_data, extra_node_ids):
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
        # write species 
        for entry in aggregated_data:
            if 'accession' in entry.keys():
                label = entry['accession']
                strategy_values = [entry[strategy] for strategy in sorted(entry.keys()) if 'Single' in strategy]
                values_str = ', '.join(map(str, strategy_values))
                f.write(f"{label}, {values_str}\n")
            else:
                label = list(entry.values())[0]
                strategy_values = [entry[strategy] for strategy in sorted(entry.keys()) if 'Single' in strategy]
                values_str = ', '.join(map(str, strategy_values))
                f.write(f"{label}, {values_str}\n")
                f.write(f"'{label}', {values_str}\n")

        for label in extra_node_ids:
            all_values = []
            parts = label.split('|')
            for part in parts:
                for entry in aggregated_data:
                    if part in entry.values():
                        strategy_values = [entry[strategy] for strategy in sorted(entry.keys()) if 'Single' in strategy]
                        values_str = ', '.join(map(str, strategy_values))
                        if len(all_values) == 0:
                            all_values = values_str
                        else:
                            for i in range(len(all_values)):
                                if values_str[i] == 1:
                                    all_values[i] = 1
            f.write(f"{label}, {all_values}\n")
            f.write(f"'{label}', {all_values}\n")

def write_multis(aggregated_data, extra_node_ids):
    itol_header = """DATASET_BINARY
SEPARATOR COMMA
DATASET_LABEL,Multiple
COLOR,#aaaaaa
FIELD_COLORS,#0072B2
FIELD_LABELS,1
FIELD_SHAPES,2
DASHED_LINES,1
DATA"""
    with open('itol_multis.txt', 'w') as f:
        f.write(itol_header + '\n')
        # write species 
        for entry in aggregated_data:
            if 'accession' in entry.keys():
                label = entry['accession']
                strategy_values = [entry[strategy] for strategy in sorted(entry.keys()) if 'Multiple' in strategy]
                values_str = ', '.join(map(str, strategy_values))
                f.write(f"{label}, {values_str}\n")
            else:
                label = list(entry.values())[0]
                strategy_values = [entry[strategy] for strategy in sorted(entry.keys()) if 'Multiple' in strategy]
                values_str = ', '.join(map(str, strategy_values))
                f.write(f"{label}, {values_str}\n")
                f.write(f"'{label}', {values_str}\n")

        for label in extra_node_ids:
            all_values = []
            parts = label.split('|')
            for part in parts:
                for entry in aggregated_data:
                    if part in entry.values():
                        strategy_values = [entry[strategy] for strategy in sorted(entry.keys()) if 'Multiple' in strategy]
                        values_str = ', '.join(map(str, strategy_values))
                        if len(all_values) == 0:
                            all_values = values_str
                        else:
                            for i in range(len(all_values)):
                                if values_str[i] == 1:
                                    all_values[i] = 1
            f.write(f"{label}, {all_values}\n")
            f.write(f"'{label}', {all_values}\n")

def write_combo12(aggregated_data, extra_node_ids):
    itol_header = """DATASET_PIECHART
SEPARATOR COMMA
DATASET_LABEL,Combo 1+2 clusters - pie
COLOR,#aaaaaa
FIELD_COLORS,#56B4E9,#D55E00,#E69F00,#F0E442
FIELD_LABELS,1,2,3,4
DASHED_LINES,0
DATA"""
    with open('itol_combo12.txt', 'w') as f:
        f.write(itol_header + '\n')
        for entry in aggregated_data:
            if 'accession' in entry.keys():
                label = entry['accession']
                if entry['Combination 1 + 2'] == 1:
                    values_str = '-1,10,1,1,0,0'
                else:
                    values_str = '-1,10,0,0,0,0'
                f.write(f"{label}, {values_str}\n")
            else:
                label = list(entry.values())[0]
                if entry['Combination 1 + 2'] == 1:
                    values_str = '-1,10,1,1,0,0'
                else:
                    values_str = '-1,10,0,0,0,0'
                f.write(f"{label}, {values_str}\n")
                f.write(f"'{label}', {values_str}\n")

        for label in extra_node_ids:
            all_values = []
            parts = label.split('|')
            for part in parts:
                for entry in aggregated_data:
                    if part in entry.values():
                        if entry['Combination 1 + 2'] == 1:
                            values_str = '-1,10,1,1,0,0'
                        else:
                            values_str = '-1,10,0,0,0,0'
                        if len(all_values) == 0:
                            all_values = values_str
                        else:
                            for i in range(len(all_values)):
                                if values_str[i] == 1:
                                    all_values[i] = 1
            f.write(f"{label}, {all_values}\n")
            f.write(f"'{label}', {all_values}\n")

def write_combo34(aggregated_data, extra_node_ids):
    itol_header = """DATASET_PIECHART
SEPARATOR COMMA
DATASET_LABEL,Combo 3+4 clusters - pie
COLOR,#aaaaaa
FIELD_COLORS,#56B4E9, #D55E00, #E69F00, #F0E442
FIELD_LABELS,1,2,3,4
DASHED_LINES,0
DATA"""
    with open('itol_combo34.txt', 'w') as f:
        f.write(itol_header + '\n')
        # write species 
        for entry in aggregated_data:
            if 'accession' in entry.keys():
                label = entry['accession']
                if entry['Combination 3 + 4'] == 1:
                    values_str = '-1,10,0,0,1,1'
                else:
                    values_str = '-1,10,0,0,0,0'
                f.write(f"{label}, {values_str}\n")
            else:
                label = list(entry.values())[0]
                if entry['Combination 3 + 4'] == 1:
                    values_str = '-1,10,0,0,1,1'
                else:
                    values_str = '-1,10,0,0,0,0'
                f.write(f"{label}, {values_str}\n")
                f.write(f"'{label}', {values_str}\n")

        for label in extra_node_ids:
            all_values = []
            parts = label.split('|')
            for part in parts:
                for entry in aggregated_data:
                    if part in entry.values():
                        if entry['Combination 3 + 4'] == 1:
                            values_str = '-1,10,0,0,1,1'
                        else:
                            values_str = '-1,10,0,0,0,0'
                        if len(all_values) == 0:
                            all_values = values_str
                        else:
                            for i in range(len(all_values)):
                                if values_str[i] == 1:
                                    all_values[i] = 1
            f.write(f"{label}, {all_values}\n")
            f.write(f"'{label}', {all_values}\n")

if __name__ == "__main__":
    main()
