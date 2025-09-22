#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Plot a violin plot of protein pI.')
    parser.add_argument('--df', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--out', type=str, default='../outputs', help='Directory to save the plot')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    df = pd.read_csv(args.df)

    if 'pI' not in df.columns:
        print("Error: 'pI' column not found in input file.")
        return

    plt.figure(figsize=(4, 3))
    sns.violinplot(y='pI', data=df, color='#999999', cut=0)
    plt.title('Histone pI')
    plt.ylim(0,14)
    plt.ylabel('pI')
    plt.tight_layout()
    plt.savefig(os.path.join(args.out, 'protein_pI_violin.png'), dpi=300)
    plt.close()

    # Output how many histones are under pI 7
    count_under_7 = (df['pI'] < 7).sum()
    count_over_7 = (df['pI'] > 7).sum()
    print(f"Number of histones with pI < 7: {count_under_7}")
    print(f"Number of histones with pI > 7: {count_over_7}")

if __name__ == "__main__":
    main()