#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
import argparse
import re


def main():
    # Set up argparse to take in the CSV file as an argument
    parser = argparse.ArgumentParser(description='Plot strategy distribution from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    args = parser.parse_args()
    plot_strategy_distribution(args.df)
    
def plot_strategy_distribution(csv_file):
    df = pd.read_csv(csv_file)
    strategy_counts = df['Strategy'].value_counts()

    # Prepare data for plotting
    clusters = strategy_counts.index.tolist()
    values = strategy_counts.values.tolist()

    # Define colors for the bars
    hex_colors = ['#aec7e8', '#ffb343', '#98df8a', '#ff9896', '#c5b0d5', '#c49c8a', 
                  '#fbb4e0', '#cccccc', '#dbdb8d', '#9edae5']
    colors = []
    
    for cluster in clusters:
        if '+' not in cluster:
            # Extract the last digit from a cluster name like 'Single 2' or 'Multiple 1'
            match = re.search(r'\d+', cluster)
            if match:
                index = int(match.group())  # Convert the extracted number to an integer
                colors.append(hex_colors[index + 1])  # Add 1 for the color index
            else:
                colors.append('#aaaaaa')  # Fallback color if no number is found
        else:
            colors.append('#aaaaaa')  # Gray for combinations

    # Create the bar plot
    fig, ax = plt.subplots()
    bars = ax.bar(clusters, values, color=colors)
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters, rotation=90)

    # Add value labels on top of bars
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, 
                 bar.get_height(), 
                 f'{value}', 
                 ha='center', va='bottom')

    # Set title and remove top/right spines
    ax.set_title('Strategy Used by Each Genome')
    ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig('../outputs/strategies.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
