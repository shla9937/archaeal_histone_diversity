#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import argparse
import ast

HYDROPHOBIC = ['A', 'I', 'L', 'V', 'M', 'P']
AROMATIC = ['Y', 'F', 'W']
CHARGED = ['K', 'R', 'E', 'D']
HISTIDINES = ['H']
ALL_AMINO_ACIDS = list('GSTQNC')

def categorize_amino_acids():
    sorted_aa = list(HISTIDINES + HYDROPHOBIC + AROMATIC + CHARGED + ALL_AMINO_ACIDS)
    return sorted_aa

def main():
    parser = argparse.ArgumentParser(description="Plot amino acid distributions from histone sequences for different strategies.")
    parser.add_argument('--strategies_df', type=str, required=True, help="Path to the CSV file with strategies and histones.")
    parser.add_argument('--sequences_df', type=str, required=True, help="Path to the CSV file with protein sequences.")
    args = parser.parse_args()

    global cud_palette
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"]

    strategies_df = pd.read_csv(args.strategies_df)
    sequences_df = pd.read_csv(args.sequences_df)
    strategies_df['Clusters'] = strategies_df['Clusters'].apply(ast.literal_eval)
    strategies_df['Histones'] = strategies_df['Histones'].apply(ast.literal_eval)
    
    for strategy in strategies_df['Strategy'].unique():
        strategy_df = strategies_df[strategies_df['Strategy'] == strategy]
        if 'Combination' in strategy:
            plot_combination_histone_distribution(strategy_df, strategy, sequences_df)
        else:
            plot_single_histone_distribution(strategy_df, strategy, sequences_df)

    plot_combined_amino_acid_distribution(strategies_df, sequences_df)

def get_sequence_for_histones(histone_list, sequences_df):
    histone_sequences = {}
    for histone in histone_list:
        seq_row = sequences_df[sequences_df['Protein Name'] == histone]
        if not seq_row.empty:
            histone_sequences[histone] = seq_row['Sequence'].values[0]
    return histone_sequences

def calculate_weighted_percentage_amino_acids(sequences):
    amino_acid_frequencies = Counter()
    total_sequences = len(sequences)
    
    if total_sequences == 0:
        return Counter()

    for seq in sequences:
        sequence_length = len(seq)
        if sequence_length > 0:
            sequence_counter = Counter(seq)
            for aa, count in sequence_counter.items():
                amino_acid_frequencies[aa] += (count / sequence_length) * 100
    return {aa: freq / total_sequences for aa, freq in amino_acid_frequencies.items()}

def plot_combination_histone_distribution(strategy_df, strategy, sequences_df):
    fig, ax = plt.subplots(figsize=(2.5, 1.5))    
    plt.rcParams.update({'font.size': 6})
    fig.subplots_adjust(bottom=0.5)
    ax.spines[['right', 'top']].set_visible(False)

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
        cluster_aa_percentages = calculate_weighted_percentage_amino_acids(sequences)
        y_values = [cluster_aa_percentages.get(aa, 0) for aa in all_aa]
        positions = [x + i * bar_width for x in x_indices]
        plt.bar(positions, y_values, width=bar_width, color=cud_palette[int(cluster)], label=f'Cluster {cluster}')

    plt.xlabel('Amino Acids', fontsize=6)
    plt.ylabel('Percentage (%)', fontsize=6)
    plt.title(f'AA Distribution for {strategy}')
    plt.xticks([x + bar_width * (len(cluster_map) - 1) / 2 for x in x_indices], all_aa, fontsize=6)
    plt.yticks(fontsize=6)
    plt.ylim(0, 25)
    plt.legend(loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'../outputs/{strategy}_amino_acid_distribution.png', dpi=300)
    plt.show()

def plot_single_histone_distribution(strategy_df, strategy, sequences_df):
    fig, ax = plt.subplots(figsize=(2.5, 1.5))    
    plt.rcParams.update({'font.size': 6})
    fig.subplots_adjust(bottom=0.5)
    ax.spines[['right', 'top']].set_visible(False)

    histone_aa_counts = Counter()
    histone_list = strategy_df['Histones'].explode().unique()
    histone_sequences = get_sequence_for_histones(histone_list, sequences_df)
    all_sequences = list(histone_sequences.values())
    histone_aa_percentages = calculate_weighted_percentage_amino_acids(all_sequences)

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

    plt.xlabel('Amino Acids', fontsize=6)
    plt.ylabel('Percentage (%)', fontsize=6)
    plt.title(f'AA Distribution for {strategy}')
    plt.xticks(rotation=0, fontsize=6)
    plt.ylim(0, 25)
    plt.yticks(fontsize=6)
    for x, y in zip(x_values, y_values):
        if y > 0:
            plt.text(x, y, f'{y:.2f}%', ha='center', va='bottom', rotation=90)
    plt.tight_layout()
    plt.savefig(f'../outputs/{strategy}_single_amino_acid_distribution.png', dpi=300)
    plt.show()

def plot_combined_amino_acid_distribution(strategies_df, sequences_df):
    fig, ax = plt.subplots(figsize=(2.5, 1.5))    
    plt.rcParams.update({'font.size': 6})
    fig.subplots_adjust(bottom=0.5)
    ax.spines[['right', 'top']].set_visible(False)

    all_sequences = []
    for _, row in strategies_df.iterrows():
        histone_list = row['Histones']
        histone_sequences = get_sequence_for_histones(histone_list, sequences_df)
        all_sequences.extend(histone_sequences.values())
    
    combined_aa_percentages = calculate_weighted_percentage_amino_acids(all_sequences)
    x_values = categorize_amino_acids()
    y_values = [combined_aa_percentages.get(aa, 0) for aa in x_values]

    plt.bar(x_values, y_values, color='#aaaaaa')  # Use a distinct color for the combined plot
    plt.xlabel('Amino Acids', fontsize=6)
    plt.ylabel('Percentage (%)', fontsize=6)
    plt.title('AA Distribution Across All Strategies')
    plt.xticks(rotation=0, fontsize=6)
    plt.ylim(0, 25)
    plt.yticks(fontsize=6)
    for x, y in zip(x_values, y_values):
        if y > 0:
            plt.text(x, y, f'{y:.2f}%', ha='center', va='bottom', rotation=90)

    plt.tight_layout()
    plt.savefig('../outputs/combined_amino_acid_distribution.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
