#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import ast
import argparse

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Plot distribution of cluster combinations within combination and multiple strategies.")
    parser.add_argument('--df', type=str, required=True, help="Path to the CSV file.")
    args = parser.parse_args()

    global cud_palette
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"]

    # Load the DataFrame
    df = pd.read_csv(args.df)

    # Filter for only 'Combination' and 'Multiple' strategies
    combo_multi_df = df[df['Strategy'].str.contains('Combination|Multiple')]

    # Determine the maximum x-value across all strategies
    max_x_value = get_max_x_value(combo_multi_df)
    
    # Iterate over unique strategies
    for strategy in combo_multi_df['Strategy'].unique():
        strategy_df = combo_multi_df[combo_multi_df['Strategy'] == strategy]
        if 'Combination' in strategy:
            plot_combination_distribution(strategy_df, strategy, max_x_value)
        elif 'Multiple' in strategy:
            plot_multiple_distribution(strategy_df, strategy, max_x_value)

def get_max_x_value(df):
    max_x = 1
    for strategy in df['Strategy'].unique():
        strategy_df = df[df['Strategy'] == strategy]
        if 'Combination' in strategy:
            clusters_list = strategy_df['Clusters'].apply(ast.literal_eval)
            unique_clusters = sorted(set(cluster for sublist in clusters_list for cluster in sublist))
            for cluster in unique_clusters:
                cluster_counts = clusters_list.apply(lambda x: x.count(cluster))
                max_x = max(max_x, max(cluster_counts, default=0))
        elif 'Multiple' in strategy:
            clusters_list = strategy_df['Clusters'].apply(ast.literal_eval)
            cluster_counts = clusters_list.apply(len)
            max_x = max(max_x, max(cluster_counts, default=0))
    return max_x

def plot_combination_distribution(strategy_df, strategy, max_x_value):
    # Extract clusters as lists of integers
    clusters_list = strategy_df['Clusters'].apply(ast.literal_eval)
    
    # Identify unique clusters in the strategy
    unique_clusters = sorted(set(cluster for sublist in clusters_list for cluster in sublist))

    for cluster in unique_clusters:
        # Count occurrences of the current cluster
        cluster_counts = clusters_list.apply(lambda x: x.count(cluster))
        counts = Counter(cluster_counts)

        # Create a bar plot
        plt.figure(figsize=(10, 6))
        x_values = sorted(counts.keys())
        y_values = [counts[x] for x in x_values]
        plt.bar(x_values, y_values, color=cud_palette[int(cluster)])
        plt.xlabel(f'Number of Cluster {cluster} Occurrences')
        plt.ylabel('Frequency')
        plt.title(f'Distribution of Cluster {cluster} for {strategy} Strategy')
        plt.xticks(range(1, max_x_value + 1))  # Set x-axis to start at 1
        plt.xlim(0.5, max_x_value + 0.5)  # Ensure x-axis starts at 1 and extends to max_x_value

        # Add counts on top of bars
        for x, y in zip(x_values, y_values):
            plt.text(x, y, f'{y}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(f'../outputs/{strategy}_cluster_{cluster}_distribution.png', dpi=300)
        plt.show()

def plot_multiple_distribution(strategy_df, strategy, max_x_value):
    # Extract clusters as lists of integers
    clusters_list = strategy_df['Clusters'].apply(ast.literal_eval)

    # Count the number of clusters for each entry
    cluster_counts = clusters_list.apply(len)
    counts = Counter(cluster_counts)

    # Determine color based on strategy name
    strategy_index = int(strategy.split(' ')[-1])  # Get the number after "Multiple"
    color = cud_palette[strategy_index] if strategy_index < len(cud_palette) else '#aaaaaa'

    # Create a bar plot
    plt.figure(figsize=(10, 6))
    x_values = list(range(1, max_x_value + 1))  # Ensure x-values start from 1 and include all integers up to the max value
    y_values = [counts.get(x, 0) for x in x_values]
    plt.bar(x_values, y_values, color=color)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of Number of Clusters for {strategy} Strategy')
    
    plt.xticks(range(1, max_x_value + 1))  # Set x-axis to start at 1
    plt.xlim(0.5, max_x_value + 0.5)  # Ensure x-axis starts at 1 and extends to max_x_value

    # Add counts on top of bars
    for x, y in zip(x_values, y_values):
        plt.text(x, y, f'{y}', ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(f'../outputs/{strategy}_number_clusters_distribution.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
