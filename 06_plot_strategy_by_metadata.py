#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from scipy.stats import ttest_ind, shapiro, mannwhitneyu
import re

def main():
    parser = argparse.ArgumentParser(description='Plot strategy distribution from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--cutoff', type=int, required=True, help='Minimum number of genomes for a strategy to be included')
    args = parser.parse_args()

    global cud_palette
    cud_palette = ["#999999","#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000"]    

    proteins = pd.read_csv(args.df)
    proteins_exploded = explode_list_column(proteins, 'Strategy')
    strategy_counts = proteins_exploded['Strategy'].value_counts()
    valid_strategies = strategy_counts[strategy_counts > args.cutoff].index.tolist()
    filtered_df = proteins_exploded[proteins_exploded['Strategy'].isin(valid_strategies)]

    category_order = get_categories(proteins)
    valid_category_order = [strategy for strategy in category_order if strategy in valid_strategies]
    filtered_df['Strategy'] = pd.Categorical(
        filtered_df['Strategy'], categories=valid_category_order, ordered=True)
    columns_to_compare = ['Phylum', 'coding_density', 'gc_percentage', 'genome_size']
    for col in columns_to_compare:
        if col == 'Phylum':
            plot_comparison(filtered_df, 'Strategy', col, None, save_path=f'../outputs/strategy_v_{col}.png', use_colors=False)
            plot_comparison(filtered_df, 'Strategy', col, None, save_path=f'../outputs/{col}_v_strategy.png', use_colors=True)
        else:
            p_values_df = perform_pairwise_comparison(filtered_df, col)
            plot_comparison(filtered_df, 'Strategy', col, p_values_df, save_path=f'../outputs/strategy_v_{col}.png', use_colors=True)

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

# Function to plot different charts
def plot_comparison(df, x_col, y_col, p_values_df, save_path=None, use_colors=True):
    plt.figure(figsize=(16, 8))  # Increase the figure size

    if use_colors:
        colors = []

        for cluster in df[x_col].cat.categories:
            if '+' not in cluster:
                match = re.search(r'\d+', cluster)
                if 'Multiple' in cluster:
                    colors.append(cud_palette[-2])
                elif match:
                    index = int(match.group())
                    colors.append(cud_palette[index])
                else:
                    colors.append('#444444')
            else:
                colors.append('#aaaaaa')
    else:
        colors = None

    if y_col == 'Phylum' and use_colors == True:
        sns.countplot(data=df, x=y_col, hue=x_col, palette=colors)
        plt.xticks(rotation=90)
    elif y_col == 'Phylum':
        sns.countplot(data=df, x=x_col, hue=y_col)
        plt.xticks(rotation=90)
    else:
        sns.boxplot(data=df, x=x_col, y=y_col, palette=colors)
        plt.xticks(rotation=90)

        # Add significant p-values to the plot
        if p_values_df is not None:
            for i, strategy in enumerate(df[x_col].cat.categories):
                if strategy != 'No histones':
                    p_value = p_values_df.loc[p_values_df['Strategy'] == strategy, 'p-value'].values[0]
                    if p_value < 0.05:
                        plt.text(i, df[y_col].max(), f'{p_value:.3e}', ha='center', color='red')

    plt.title(f'{x_col} vs {y_col}')
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.tight_layout()
    
    # Adjust the legend to fit better in the plot area
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title=y_col, fontsize='small', title_fontsize='medium')  # Adjust legend size
    plt.subplots_adjust(right=0.85)  # Adjust the right space to fit the legend

    if save_path:
        save_path = save_path.replace(' ', '_')  # Remove spaces from file name
        plt.savefig(save_path, bbox_inches='tight')  # Use bbox_inches to avoid cutoff  
    plt.show()

def perform_pairwise_comparison(df, y_col):
    # Extract unique strategies
    strategies = df['Strategy'].unique()

    # Create an empty DataFrame to store the values
    comparison_df = pd.DataFrame()

    for strategy in strategies:
        strategy_data = df[df['Strategy'] == strategy][y_col].reset_index(drop=True)
        comparison_df[strategy] = strategy_data

    # Fill NaNs with None for consistent output
    comparison_df = comparison_df.where(pd.notnull(comparison_df), None)

    # Perform pairwise tests
    results = []
    for strategy in strategies:
        if strategy != 'No histones':
            strategy_data = comparison_df[strategy].dropna().astype(float)
            no_histones_data = comparison_df['No histones'].dropna().astype(float)
            
            if not strategy_data.empty and not no_histones_data.empty:
                _, p_no_histones = shapiro(no_histones_data)
                _, p_strategy = shapiro(strategy_data)
                
                if p_no_histones > 0.05 and p_strategy > 0.05:  # If both are normally distributed
                    stat, p_value = ttest_ind(no_histones_data, strategy_data, equal_var=False)  # Welch's t-test
                    results.append({
                        'Strategy': strategy,
                        't-statistic': stat,
                        'p-value': p_value
                    })
                else:
                    stat, p_value = mannwhitneyu(no_histones_data, strategy_data, alternative='two-sided')  # Mann-Whitney U test
                    results.append({
                        'Strategy': strategy,
                        't-statistic': stat,
                        'p-value': p_value
                    })
            else:
                results.append({
                    'Strategy': strategy,
                    't-statistic': None,
                    'p-value': None
                })

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    return results_df

if __name__ == "__main__":
    main()
