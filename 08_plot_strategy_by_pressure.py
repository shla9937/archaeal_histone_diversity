#!/usr/bin/env python

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

    global cud_palette
    cud_palette = ["#999999", "#56B4E9", "#D55E00", "#E69F00", "#F0E442", "#009E73", "#CC79A7", "#0072B2", "#000000"]

    species_location = pd.read_csv(args.df)
    species_location_exploded = explode_list_column(species_location, 'Strategy')
    strategy_counts = species_location_exploded['Strategy'].value_counts()
    valid_strategies = strategy_counts[strategy_counts > args.cutoff].index.tolist()
    filtered_df = species_location_exploded[species_location_exploded['Strategy'].isin(valid_strategies)]
    
    # Generate individual proportional strategy plots
    plot_pressure_per_strategy(filtered_df, save_dir='../outputs/')
    
    # Generate combined proportion plot
    plot_combined_pressure_proportions(filtered_df, save_dir='../outputs/')

def explode_list_column(df, col):
    return df.explode(col)

# Function to assign color based on strategy
def get_color_for_strategy(strategy):
    if 'Single' in strategy:
        cluster_number = strategy.replace('Single', '').strip()
        return cud_palette[int(cluster_number) % len(cud_palette)]
    elif 'Multiple' in strategy:
        cluster_number = strategy.replace('Multiple', '').strip()
        return cud_palette[(-int(cluster_number) - 1) % len(cud_palette)]
    elif 'Combination' in strategy:
        return '#aaaaaa'
    elif 'No histones' in strategy:
        return cud_palette[-1]
    else:
        return cud_palette[-1]  # Default color if strategy doesn't match any case

def plot_pressure_per_strategy(df, save_dir='../outputs/'):
    os.makedirs(save_dir, exist_ok=True)
    
    pressure_data = []
    total_genomes_per_strategy = df.groupby('Strategy').size().reset_index(name='Total_genomes')
    
    for _, row in df.iterrows():
        pressures = row['pressure']
        strategy = row['Strategy']
        if isinstance(pressures, str):
            pressures = eval(pressures)  # Assuming the pressure is a list in string format
        
        # Filter out 'none' and 'other' pressures
        pressures = [p for p in pressures if p.lower() not in ['none', 'other']]
        
        for pressure in pressures:
            pressure_data.append({'Strategy': strategy, 'Pressure': pressure})
    
    pressure_df = pd.DataFrame(pressure_data)
    pressure_counts = pressure_df.groupby(['Strategy', 'Pressure']).size().reset_index(name='Count')
    
    # Merge with total genomes to calculate proportions
    pressure_proportions = pressure_counts.merge(total_genomes_per_strategy, on='Strategy')
    pressure_proportions['Proportion'] = (pressure_proportions['Count'] / pressure_proportions['Total_genomes']) * 100  # Convert to percentage
    plt.rcParams.update({'font.size': 6})
    for strategy in pressure_proportions['Strategy'].unique():
        strategy_data = pressure_proportions[pressure_proportions['Strategy'] == strategy]
        strategy_data = strategy_data.sort_values(by='Proportion', ascending=False)

        # Ensure consistent figure size and formatting each loop iteration
        fig, ax = plt.subplots(figsize=(2.5, 2))
        fig.subplots_adjust(bottom=0.4)   
        # Get the color for the current strategy
        color = get_color_for_strategy(strategy)
        
        # Plot with consistent color palette and formatting
        sns.barplot(data=strategy_data, x='Pressure', y='Proportion', palette=[color]*len(strategy_data), ax=ax)
        ax.spines[['right', 'top']].set_visible(False)
        ax.set_xlabel('Pressure', fontsize=6)
        plt.rcParams.update({'font.size': 6})
        ax.set_ylabel('Percentage of genomes', fontsize=6)        
        plt.title(f'Pressures for {strategy}')
        plt.xticks(rotation=90)
        plt.ylim(0, 100)  # Set y-axis limits to 0-100%

        # Tight layout ensures no clipping
        plt.tight_layout()
        
        # Save the figure and then close it to avoid interference with subsequent plots
        save_path = os.path.join(save_dir, f'pressure_percentage_{strategy}.png')
        plt.savefig(save_path, dpi=300)

def plot_combined_pressure_proportions(df, save_dir='../outputs/'):
    os.makedirs(save_dir, exist_ok=True)
    
    # Aggregate data for proportions
    pressure_data = []
    total_genomes_per_strategy = df.groupby('Strategy').size().reset_index(name='Total_genomes')
    
    for _, row in df.iterrows():
        pressures = row['pressure']
        strategy = row['Strategy']
        if isinstance(pressures, str):
            pressures = eval(pressures)  # Assuming the pressure is a list in string format
        
        # Filter out 'none' and 'other' pressures
        pressures = [p for p in pressures if p.lower() not in ['none', 'other']]
        
        for pressure in pressures:
            pressure_data.append({'Strategy': strategy, 'Pressure': pressure})
    
    pressure_df = pd.DataFrame(pressure_data)
    pressure_counts = pressure_df.groupby(['Strategy', 'Pressure']).size().reset_index(name='Count')
    
    # Merge with total genomes to calculate proportions
    pressure_proportions = pressure_counts.merge(total_genomes_per_strategy, on='Strategy')
    pressure_proportions['Proportion'] = (pressure_proportions['Count'] / pressure_proportions['Total_genomes']) * 100  # Convert to percentage
    
    # Define the desired order for strategies: Singles, Multiples, Combinations, No histones
    strategy_order = (
        [s for s in pressure_proportions['Strategy'].unique() if 'Single' in s] +
        [s for s in pressure_proportions['Strategy'].unique() if 'Multiple' in s] +
        [s for s in pressure_proportions['Strategy'].unique() if 'Combination' in s] +
        [s for s in pressure_proportions['Strategy'].unique() if 'No histones' in s]
    )
    
    # Sort the data based on this order
    pressure_proportions['Strategy'] = pd.Categorical(pressure_proportions['Strategy'], categories=strategy_order, ordered=True)
    pressure_proportions = pressure_proportions.sort_values(['Strategy', 'Pressure'])
    
    # Create a color palette based on the pressure
    unique_pressures = pressure_proportions['Pressure'].unique()
    
    # Plot proportions for each strategy as x and group the bars by pressure
    fig, ax = plt.subplots(figsize=(2.5, 2))    
    plt.rcParams.update({'font.size': 6})
    fig.subplots_adjust(bottom=0.5)
    
    ax = sns.barplot(data=pressure_proportions, x='Strategy', y='Proportion', hue='Pressure', 
                     hue_order=unique_pressures,
                     palette='tab10')  # Adjust color palette as needed
    
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel('Strategy', fontsize=6)
    ax.set_ylabel('Percentage of genomes', fontsize=6)

    plt.title('Proportion of Pressures by Strategy')
    plt.xticks(rotation=90)
    plt.ylim(0, 100)  # Set y-axis limits to 0-100%
    ax.legend_.set_title('Pressure')
    plt.tight_layout()
    
    save_path = os.path.join(save_dir, 'combined_pressure_percentages_by_strategy.png')
    plt.savefig(save_path, dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
    