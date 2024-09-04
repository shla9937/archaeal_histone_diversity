#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Plot strategy vs pressure from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--cutoff', type=int, required=True, help='Minimum number of genomes for a strategy to be included')
    args = parser.parse_args()

    species_location = pd.read_csv(args.df)
    species_location_exploded = explode_list_column(species_location, 'Strategy')
    strategy_counts = species_location_exploded['Strategy'].value_counts()
    valid_strategies = strategy_counts[strategy_counts > args.cutoff].index.tolist()
    filtered_df = species_location_exploded[species_location_exploded['Strategy'].isin(valid_strategies)]
    plot_pressure_per_strategy(filtered_df, save_dir='../outputs/')

def explode_list_column(df, col):
    return df.explode(col)

def plot_pressure_per_strategy(df, save_dir='../outputs/'):
    os.makedirs(save_dir, exist_ok=True)
    pressure_data = []
    for _, row in df.iterrows():
        pressures = row['pressure']
        strategy = row['Strategy']
        if isinstance(pressures, str):
            pressures = eval(pressures)
        for pressure in pressures:
            pressure_data.append({'Strategy': strategy, 'Pressure': pressure})
    pressure_df = pd.DataFrame(pressure_data)
    pressure_counts = pressure_df.groupby(['Strategy', 'Pressure']).size().reset_index(name='Count')
    for strategy in pressure_counts['Strategy'].unique():
        strategy_data = pressure_counts[pressure_counts['Strategy'] == strategy]
        strategy_data = strategy_data.sort_values(by='Count', ascending=False)
        plt.figure(figsize=(10, 6))
        ax = sns.barplot(data=strategy_data, x='Pressure', y='Count', palette='viridis')
        plt.title(f'Pressure Distribution for {strategy}')
        plt.xlabel('Pressure')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout()
        save_path = os.path.join(save_dir, f'pressure_distribution_{strategy}.png')
        plt.savefig(save_path)
        plt.show()

if __name__ == "__main__":
    main()
