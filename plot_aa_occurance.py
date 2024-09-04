#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import argparse
import ast

# Define amino acid categories
HYDROPHOBIC = ['A', 'I', 'L', 'V', 'M', 'P']
AROMATIC = ['Y', 'F', 'W']
CHARGED = ['K', 'R', 'E', 'D']
HISTIDINES = ['H']
ALL_AMINO_ACIDS = list('GSTQNC')

def categorize_amino_acids():
    """Categorize amino acids."""
    sorted_aa = list(
        HISTIDINES +
        HYDROPHOBIC +
        AROMATIC +
        CHARGED +
        ALL_AMINO_ACIDS
    )
    return sorted_aa

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Plot amino acid distributions from histone sequences for different strategies.")
    parser.add_argument('--df1', type=str, required=True, help="Path to the CSV file with strategies and histones.")
    parser.add_argument('--df2', type=str, required=True, help="Path to the CSV file with protein sequences.")
    args = parser.parse_args()

    global cud_palette
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"]

    # Load the DataFrames
    strategies_df = pd.read_csv(args.df1)
    sequences_df = pd.read_csv(args.df2)

    # Convert 'Clusters' and 'Histones' columns to lists
    strategies_df['Clusters'] = strategies_df['Clusters'].apply(ast.literal_eval)
    strategies_df['Histones'] = strategies_df['Histones'].apply(ast.literal_eval)

    # Process each strategy
    for strategy in strategies_df['Strategy'].unique():
        strategy_df = strategies_df[strategies_df['Strategy'] == strategy]

        # For combinations, plot side-by-side bars for histone types
        if 'Combination' in strategy:
            plot_combination_histone_distribution(strategy_df, strategy, sequences_df)
        else:
            plot_single_histone_distribution(strategy_df, strategy, sequences_df)

def get_sequence_for_histones(histone_list, sequences_df):
    """Retrieve sequences for given histones."""
    histone_sequences = {}
    for histone in histone_list:
        seq_row = sequences_df[sequences_df['Protein Name'] == histone]
        if not seq_row.empty:
            histone_sequences[histone] = seq_row['Sequence'].values[0]
    return histone_sequences

def calculate_percentage_amino_acids(sequences):
    """Calculate the percentage of each amino acid in the sequences."""
    total_aa = sum(len(seq) for seq in sequences)
    if total_aa == 0:
        return Counter()
    amino_acid_counter = Counter()
    for seq in sequences:
        amino_acid_counter.update(seq)
    return {aa: (count / total_aa) * 100 for aa, count in amino_acid_counter.items()}

def plot_combination_histone_distribution(strategy_df, strategy, sequences_df):
    plt.figure(figsize=(14, 8))

    # Determine unique clusters from the combinations
    cluster_map = {}
    for clusters, histones in zip(strategy_df['Clusters'], strategy_df['Histones']):
        for cluster, histone_list in zip(clusters, histones):
            if cluster not in cluster_map:
                cluster_map[cluster] = []
            histone_sequences = get_sequence_for_histones([histone_list], sequences_df)
            cluster_map[cluster].extend(histone_sequences.values())

    all_aa = categorize_amino_acids()
    bar_width = 0.35
    x_indices = range(len(all_aa))

    for i, (cluster, sequences) in enumerate(cluster_map.items()):
        cluster_aa_percentages = calculate_percentage_amino_acids(sequences)
        y_values = [cluster_aa_percentages.get(aa, 0) for aa in all_aa]
        positions = [x + i * bar_width for x in x_indices]
        plt.bar(positions, y_values, width=bar_width, color=cud_palette[int(cluster)], label=f'Cluster {cluster}')

    plt.xlabel('Amino Acids')
    plt.ylabel('Percentage (%)')
    plt.title(f'Amino Acid Distribution for {strategy} Strategy')
    plt.xticks([x + bar_width * (len(cluster_map) - 1) / 2 for x in x_indices], all_aa)
    plt.ylim(0, 25)  # Set y-axis limit from 0 to 25%
    plt.legend(loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

    plt.tight_layout()
    plt.savefig(f'../outputs/{strategy}_amino_acid_distribution.png', dpi=300)
    plt.show()

def plot_single_histone_distribution(strategy_df, strategy, sequences_df):
    plt.figure(figsize=(12, 6))

    histone_aa_counts = Counter()
    histone_list = strategy_df['Histones'].explode().unique()
    histone_sequences = get_sequence_for_histones(histone_list, sequences_df)
    all_sequences = list(histone_sequences.values())
    histone_aa_percentages = calculate_percentage_amino_acids(all_sequences)

    # Use the order defined at the top
    x_values = categorize_amino_acids()
    y_values = [histone_aa_percentages.get(aa, 0) for aa in x_values]
    if 'Single' in strategy:
        cluster_number = strategy.replace('Single', '').strip()
    elif 'Multiple' in strategy:
        cluster_number = strategy.replace('Multiple', '').strip()
        cluster_number = -int(cluster_number)-1
    elif 'No histones' in strategy:
        cluster_number = 0
    plt.bar(x_values, y_values, color=cud_palette[int(cluster_number)])

    plt.xlabel('Amino Acids')
    plt.ylabel('Percentage (%)')
    plt.title(f'Amino Acid Distribution for {strategy} Strategy')
    plt.xticks(rotation=0)  # Ensure x-axis labels are not rotated
    plt.ylim(0, 25)  # Set y-axis limit from 0 to 25%

    # Add percentages on top of bars, but only if the percentage is not zero
    for x, y in zip(x_values, y_values):
        if y > 0:
            plt.text(x, y, f'{y:.2f}%', ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(f'../outputs/{strategy}_single_amino_acid_distribution.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
