#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def main():
    # Set up argparse to take in the CSV file
    parser = argparse.ArgumentParser(description='Plot pressures from a CSV file.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    args = parser.parse_args()

    # Read the DataFrame from the CSV file
    species_location = pd.read_csv(args.df)

    # Plot the pressure distribution for all organisms
    plot_pressure_distribution(species_location, save_dir='../outputs/')

def plot_pressure_distribution(df, save_dir='../outputs/'):
    # Ensure the save directory exists
    os.makedirs(save_dir, exist_ok=True)

    # Create a DataFrame to hold the pressure data
    pressure_data = []

    # Iterate through each row and extract pressures
    for _, row in df.iterrows():
        pressures = row['pressure']
        
        # Convert pressures string to list (if needed)
        if isinstance(pressures, str):
            pressures = eval(pressures)  # Use eval carefully if you're sure of the input format
        
        for pressure in pressures:
            pressure_data.append({'Pressure': pressure})

    # Convert the list of dictionaries to a DataFrame
    pressure_df = pd.DataFrame(pressure_data)

    # Count occurrences of each pressure type
    pressure_counts = pressure_df['Pressure'].value_counts().reset_index()
    pressure_counts.columns = ['Pressure', 'Count']

    # Sort the data by count in descending order
    pressure_counts = pressure_counts.sort_values(by='Count', ascending=False)

    # Plotting
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(data=pressure_counts, x='Pressure', y='Count', palette='viridis')
    plt.title('Pressure Distribution for All Organisms')
    plt.xlabel('Pressure')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    
    # Save the plot to the outputs directory
    save_path = os.path.join(save_dir, 'pressure_distribution_all_organisms.png')
    plt.savefig(save_path)
    plt.show()

if __name__ == "__main__":
    main()
