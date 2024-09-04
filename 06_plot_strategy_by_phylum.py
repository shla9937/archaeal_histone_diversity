#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from scipy.stats import ttest_ind, shapiro, mannwhitneyu
import re

def main():
    parser = argparse.ArgumentParser(description='Plot vertical boxplots for strategy comparison from a CSV file.')
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

    columns_to_compare = ['coding_density', 'gc_percentage', 'genome_size']
    for col in columns_to_compare:
        p_values_df = perform_pairwise_comparison(filtered_df, col)
        save_p_values(p_values_df, col)
        plot_comparison(filtered_df, 'Strategy', col, p_values_df, save_path=f'../outputs/strategy_v_{col}_boxplot.png', use_colors=True)

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
    
    plt.figure(figsize=(10, 8))  # Adjusted the figure size for better fit

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

    sns.boxplot(data=df, x=x_col, y=y_col, palette=colors, orient='v')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability

    # Set y-axis limits for percentage columns
    if y_col in ['gc_percentage', 'coding_density']:
        plt.ylim(0, 105)
    else:
        plt.ylim(df[y_col].min() * 0.9, df[y_col].max() * 1.1)  # Default scaling

    # Add significance stars at a fixed height
    if p_values_df is not None:
        strategies = df[x_col].cat.categories
        no_histones_index = strategies.get_loc('No histones')
        y_max = plt.ylim()[1]  # Set to the maximum y limit
        
        # Choose a fixed height for stars
        star_height = y_max - 0.05 * y_max  # 5% below the maximum y limit

        # Adding stars
        for i, strategy in enumerate(strategies):
            if strategy != 'No histones':
                p_value = p_values_df.loc[p_values_df['Strategy'] == strategy, 'p-value'].values[0]
                if p_value < 0.001:
                    star = '***'
                elif p_value < 0.01:
                    star = '**'
                elif p_value < 0.05:
                    star = '*'
                else:
                    star = ''
                
                plt.text(i, star_height, star, ha='center', color='red', fontsize=12)

    plt.title(f'{x_col} vs {y_col} (Comparison to No Histones)')
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.tight_layout()
    
    if save_path:
        save_path = save_path.replace(' ', '_')
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


def perform_pairwise_comparison(df, y_col):
    strategies = df['Strategy'].unique()
    comparison_df = pd.DataFrame()

    for strategy in strategies:
        strategy_data = df[df['Strategy'] == strategy][y_col].reset_index(drop=True)
        comparison_df[strategy] = strategy_data
    comparison_df = comparison_df.where(pd.notnull(comparison_df), None)
    results = []
    no_histones_data = comparison_df['No histones'].dropna().astype(float)

    for strategy in strategies:
        if strategy != 'No histones':
            strategy_data = comparison_df[strategy].dropna().astype(float)
            
            if not strategy_data.empty and not no_histones_data.empty:
                _, p_no_histones = shapiro(no_histones_data)
                _, p_strategy = shapiro(strategy_data)
                
                if p_no_histones > 0.05 and p_strategy > 0.05:
                    stat, p_value = ttest_ind(no_histones_data, strategy_data, equal_var=False)
                    results.append({
                        'Strategy': strategy,
                        't-statistic': stat,
                        'p-value': p_value
                    })
                else:
                    stat, p_value = mannwhitneyu(no_histones_data, strategy_data, alternative='two-sided')
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
    results_df = pd.DataFrame(results)
    return results_df

def save_p_values(p_values_df, y_col):
    output_file = f'../outputs/p_values_{y_col}.txt'
    with open(output_file, 'w') as f:
        f.write('Strategy\tp-value\n')
        for _, row in p_values_df.iterrows():
            f.write(f"{row['Strategy']}\t{row['p-value']:.3e}\n")

if __name__ == "__main__":
    main()
