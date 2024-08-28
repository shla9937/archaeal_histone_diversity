#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re


def main():
    parser = argparse.ArgumentParser(description='Plot strategy distribution from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    args = parser.parse_args()
    
    global cud_palette
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"]   
    
    plot_strategy_distribution(args.df)

def plot_strategy_distribution(csv_file):
    df = pd.read_csv(csv_file)
    strategy_counts = df['Strategy'].value_counts()
    clusters = strategy_counts.index.tolist()
    values = strategy_counts.values.tolist()
    
    colors = []
    for cluster in clusters:
        if 'No histones' in cluster:
            colors.append('#444444')
        elif '+' not in cluster:
            match = re.search(r'\d+', cluster)
            if 'Multiple' in cluster:
                colors.append(cud_palette[-2])
            elif match:
                index = int(match.group())
                colors.append(cud_palette[index])
            else:
                colors.append('#aaaaaa')
        else:
            colors.append('#aaaaaa')
    fig, ax = plt.subplots()
    bars = ax.bar(clusters, values, color=colors)
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters, rotation=90)
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, 
                 bar.get_height(), 
                 f'{value}', 
                 ha='center', va='bottom')
    ax.set_title('Strategies used by each genome')
    ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig('../outputs/strategies.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
