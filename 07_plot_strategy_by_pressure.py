#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def main():
    # Set up argparse to take in the CSV file and cutoff number as arguments
    parser = argparse.ArgumentParser(description='Plot strategy vs pressure from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--cutoff', type=int, required=True, help='Minimum number of genomes for a strategy to be included')
    args = parser.parse_args()

    # Read the DataFrame from the CSV file
    species_location = pd.read_csv(args.df)

    # Exploding the 'Strategy' column
    species_location_exploded = explode_list_column(species_location, 'Strategy')

    # Filter strategies based on the specified cutoff
    strategy_counts = species_location_exploded['Strategy'].value_counts()
    valid_strategies = strategy_counts[strategy_counts > args.cutoff].index.tolist()
    
    # Filter the DataFrame to only include valid strategies
    filtered_df = species_location_exploded[species_location_exploded['Strategy'].isin(valid_strategies)]

    # Plot the pressure distribution per strategy
    plot_pressure_per_strategy(filtered_df, save_dir='../outputs/')

# Function to explode the list column
def explode_list_column(df, col):
    return df.explode(col)

def plot_pressure_per_strategy(df, save_dir='../outputs/'):
    # Ensure the save directory exists
    os.makedirs(save_dir, exist_ok=True)

    # Create a DataFrame to hold the pressure data for each strategy
    pressure_data = []

    # Iterate through each row and extract pressures for each strategy
    for _, row in df.iterrows():
        pressures = row['pressure']
        strategy = row['Strategy']
        
        # Convert pressures string to list (if needed)
        if isinstance(pressures, str):
            pressures = eval(pressures)  # Use eval carefully if you're sure of the input format
        
        for pressure in pressures:
            pressure_data.append({'Strategy': strategy, 'Pressure': pressure})

    # Convert the list of dictionaries to a DataFrame
    pressure_df = pd.DataFrame(pressure_data)

    # Count occurrences of each pressure type for each strategy
    pressure_counts = pressure_df.groupby(['Strategy', 'Pressure']).size().reset_index(name='Count')

    # Plotting for each strategy
    for strategy in pressure_counts['Strategy'].unique():
        strategy_data = pressure_counts[pressure_counts['Strategy'] == strategy]

        # Sort the data by count in descending order
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
        
        # Save each plot to the outputs directory
        save_path = os.path.join(save_dir, f'pressure_distribution_{strategy}.png')
        plt.savefig(save_path)
        plt.show()  # Close the figure to avoid displaying it inline

if __name__ == "__main__":
    main()
