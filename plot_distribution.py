#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Plot the number of histones in each genome from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    args = parser.parse_args()
    
    plot_histone_distribution(args.df)

def plot_histone_distribution(csv_file):
    df = pd.read_csv(csv_file)
    
    # Ensure 'Number of histones' column is numeric
    df['Number of histones'] = pd.to_numeric(df['Number of histones'], errors='coerce')
    
    # Drop rows where 'Number of histones' is NaN
    df = df.dropna(subset=['Number of histones'])
    
    num_histones = df['Number of histones'].astype(int).tolist()
    total_genomes = len(df)
    
    if not num_histones:
        print("No valid histone data found.")
        return

    # Calculate the range for ticks
    min_value = min(num_histones)
    max_value = max(num_histones)
    
    # Ensure ticks are set for all possible values
    all_values = list(range(min_value, max_value + 1))
    values, counts = np.unique(num_histones, return_counts=True)
    counts_dict = dict(zip(values, counts))
    
    # Fill missing values with 0 count
    all_counts = [counts_dict.get(val, 0) for val in all_values]
    
    # Percentage of genomes with more than 0 histones
    num_genomes_with_histones = len([x for x in num_histones if x > 0])
    percentage_with_histones = (num_genomes_with_histones / total_genomes) * 100
    
    # Plotting
    fig, ax = plt.subplots()
    bars = ax.bar(all_values, all_counts, color='#aaaaaa')
    ax.set_xticks(all_values)
    
    for bar, count in zip(bars, all_counts):
        plt.text(bar.get_x() + bar.get_width() / 2, 
                 bar.get_height(), 
                 f'{count}', 
                 ha='center', va='bottom')
    
    ax.set_xlabel('Number of Histones')
    ax.set_ylabel('Number of Genomes')
    ax.set_title(f'Number of Histones in Each Genome \n({num_genomes_with_histones}/{total_genomes} genomes = '
                 f'{round(percentage_with_histones, 1)}%)')
    
    ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig('../outputs/histones.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
