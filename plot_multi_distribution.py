#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import ast
import argparse

def main():
    parser = argparse.ArgumentParser(description="Plot distribution of cluster combinations within combination and multiple strategies.")
    parser.add_argument('--df', type=str, required=True, help="Path to the CSV file.")
    args = parser.parse_args()

    global cud_palette
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"]

    df = pd.read_csv(args.df)
    combo_multi_df = df[df['Strategy'].str.contains('Combination|Multiple')]
    max_x_value = get_max_x_value(combo_multi_df)
    
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
    fig, ax = plt.subplots(figsize=(2, 1.5))    
    clusters_list = strategy_df['Clusters'].apply(ast.literal_eval)
    unique_clusters = sorted(set(cluster for sublist in clusters_list for cluster in sublist))
    cluster_colors = {cluster: cud_palette[int(cluster) % len(cud_palette)] for cluster in unique_clusters}

    bar_width = 0.35
    clusters_data = {cluster: Counter(clusters_list.apply(lambda x: x.count(cluster))) for cluster in unique_clusters}
    x_values = list(range(1, max_x_value + 1))

    for i, cluster in enumerate(unique_clusters):
        counts = clusters_data[cluster]
        y_values = [counts.get(x, 0) for x in x_values]
        positions = [x + i * bar_width for x in x_values]
        plt.bar(positions, y_values, width=bar_width, color=cluster_colors[cluster], label=f'Cluster {cluster}', align='center')

    plt.rcParams.update({'font.size': 6})
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel('Number of Clusters', fontsize=6)
    ax.set_ylabel('Frequency', fontsize=6)
    plt.title(f'{strategy}')
    plt.xticks([x + bar_width * (len(unique_clusters) - 1) / 2 for x in x_values], x_values)
    plt.xlim(0, max_x_value + bar_width * len(unique_clusters))
    plt.legend(loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

    for i, cluster in enumerate(unique_clusters):
        counts = clusters_data[cluster]
        y_values = [counts.get(x, 0) for x in x_values]
        positions = [x + i * bar_width for x in x_values]
        for pos, y in zip(positions, y_values):
            if y > 0:  # Only label if y > 0
                plt.text(pos, y, f'{y}', ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(f'../outputs/{strategy}_cluster_distribution.png', dpi=300)
    plt.show()

def plot_multiple_distribution(strategy_df, strategy, max_x_value):
    clusters_list = strategy_df['Clusters'].apply(ast.literal_eval)
    cluster_counts = clusters_list.apply(len)
    counts = Counter(cluster_counts)

    strategy_index = int(strategy.split(' ')[-1])
    strategy_index = -int(strategy_index)-1
    color = cud_palette[strategy_index] if strategy_index < len(cud_palette) else '#aaaaaa'

    fig, ax = plt.subplots(figsize=(2, 1.5))    
    plt.rcParams.update({'font.size': 6})
    fig.subplots_adjust(bottom=0.5)

    x_values = list(range(1, max_x_value + 1))
    y_values = [counts.get(x, 0) for x in x_values]
    plt.bar(x_values, y_values, color=color)
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel('Number of Clusters', fontsize=6)
    ax.set_ylabel('Frequency', fontsize=6)
    plt.title(f'{strategy}')

    plt.xticks(range(1, max_x_value + 1), fontsize=6)
    plt.xlim(0.5, max_x_value + 0.5)
    plt.yticks(fontsize=6)

    for x, y in zip(x_values, y_values):
        if y > 0:
            plt.text(x, y, f'{y}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(f'../outputs/{strategy}_number_clusters_distribution.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
