#!/usr/bin/env python3

import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

def main():
    # Set up argparse to take in multiple CSV files
    parser = argparse.ArgumentParser(description='Create a Venn diagram from protein CSV files.')
    parser.add_argument('--csvs', type=str, required=True, nargs='+', help='Paths to the input CSV files')
    parser.add_argument('--names', type=str, required=True, nargs='+', help='Names of datasets')
    parser.add_argument('--col', type=str, required=True, nargs='+', help='Column to venn')
    args = parser.parse_args()

    # Read the CSV files into DataFrames
    dataframes = [pd.read_csv(csv_file) for csv_file in args.csvs]
    names = [i for i in args.names]

    # Ensure there are at least two DataFrames
    if len(dataframes) < 2:
        raise ValueError('At least two CSV files are required to create a Venn diagram.')

    # Extract the sets of 'accession' values
    sets = [set(df['Protein Name']) for df in dataframes]
    
    # Create Venn diagram based on the number of sets
    plt.figure(figsize=(8, 8))
    if len(sets) == 2:
        venn2(sets, set_labels=[f'{names[i]}' for i in range(len(sets))])
    elif len(sets) == 3:
        venn3(sets, set_labels=[f'{names[i]}' for i in range(len(sets))])
    else:
        raise ValueError('This script currently supports up to 3 CSV files for Venn diagrams.')

    # Set title and show plot
    plt.title('Venn Diagram of Protein Accessions')
    plt.show()

if __name__ == "__main__":
    main()