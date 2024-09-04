#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import re

def main():
    global cud_palette
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"]

    parser = argparse.ArgumentParser(description='Plot vertical countplots for strategy comparison from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--cutoff', type=int, required=True, help='Minimum number of genomes for a strategy to be included')
    args = parser.parse_args()

    proteins = pd.read_csv(args.df)
    proteins_exploded = explode_list_column(proteins, 'Strategy')
    strategy_counts = proteins_exploded['Strategy'].value_counts()
    valid_strategies = strategy_counts[strategy_counts > args.cutoff].index.tolist()
    filtered_df = proteins_exploded[proteins_exploded['Strategy'].isin(valid_strategies)]

    category_order = get_categories(proteins)
    valid_category_order = [strategy for strategy in category_order if strategy in valid_strategies]
    filtered_df['Strategy'] = pd.Categorical(
        filtered_df['Strategy'], categories=valid_category_order, ordered=True)

    col = 'Phylum'
    phylum_counts = filtered_df['Phylum'].value_counts()
    phyla_above_cutoff = phylum_counts[phylum_counts >= args.cutoff].index.tolist()
    phyla_below_cutoff = phylum_counts[phylum_counts < args.cutoff].index.tolist()

    filtered_above_cutoff = filtered_df[filtered_df['Phylum'].isin(phyla_above_cutoff)]
    filtered_below_cutoff = filtered_df[filtered_df['Phylum'].isin(phyla_below_cutoff)]

    plot_comparison(filtered_above_cutoff, 'Strategy', col, None, save_path=f'../outputs/{col}_v_strategy_above_cutoff_countplot.png', use_colors=True)
    plot_comparison(filtered_below_cutoff, 'Strategy', col, None, save_path=f'../outputs/{col}_v_strategy_below_cutoff_countplot.png', use_colors=True)

def explode_list_column(df, col):
    return df.explode(col)

def get_categories(proteins):
    strategies = proteins['Strategy'].unique()
    singles = sorted([s for s in strategies if s.startswith('Single')])
    multiples = sorted([s for s in strategies if s.startswith('Multiple')])
    combinations = sorted([s for s in strategies if s.startswith('Combination')])
    no_histones = [s for s in strategies if s == 'No histones']
    category_order = singles + multiples + combinations + no_histones
    return category_order

def plot_comparison(df, x_col, y_col, p_values_df, save_path=None, use_colors=True, subset=None):
    if subset is not None:
        df = df[df[x_col].isin(subset)]
    
    # Sort y values by the number of points in descending order
    y_order = df[y_col].value_counts().index

    plt.figure(figsize=(8, 8))

    colors = None
    if use_colors:
        colors = []
        for cluster in df[x_col].cat.categories:
            if '+' not in cluster:
                match = re.search(r'\d+', cluster)
                if 'Multiple' in cluster:
                    colors.append("#0072B2")
                elif match:
                    index = int(match.group())
                    colors.append(cud_palette[index])
                else:
                    colors.append('#444444')
            else:
                colors.append('#aaaaaa')

    sns.countplot(data=df, y=y_col, hue=x_col, palette=colors, order=y_order)
    
    # Add separator lines with vertical space between y values
    for i in range(len(y_order) - 1):
        plt.axhline(i + 0.6, color='gray', linestyle='--', linewidth=1)
        plt.text(-0.02, i + 0.75, '', transform=plt.gca().transAxes)  # Adds extra space above each line
        plt.text(-0.02, i + 0.45, '', transform=plt.gca().transAxes)  # Adds extra space below each line

    plt.xticks(rotation=0)

    # Ensure counts are displayed as integers
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))

    plt.title(f'{x_col} vs {y_col}')
    plt.ylabel(y_col)
    plt.xlabel('Count')
    plt.tight_layout()
    
    if save_path:
        save_path = save_path.replace(' ', '_')
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
