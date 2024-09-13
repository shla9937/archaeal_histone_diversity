#!/usr/bin/env python3

import pandas as pd
import argparse
import ast

def main():
    parser = argparse.ArgumentParser(description="Add sequences to strategies_df based on histone IDs.")
    parser.add_argument('--strategies_df', type=str, required=True, help="Path to the CSV file with strategies and histones.")
    parser.add_argument('--sequences_df', type=str, required=True, help="Path to the CSV file with protein sequences.")
    args = parser.parse_args()


    strategies_df = pd.read_csv(args.strategies_df)
    sequences_df = pd.read_csv(args.sequences_df)
    strategies_df['Histones'] = strategies_df['Histones'].apply(ast.literal_eval)
    strategies_df = add_sequences_column(strategies_df, sequences_df)
    output_strategies_df(strategies_df)

def get_sequence_for_histones(histone_list, sequences_df):
    sequences = []
    for histone in histone_list:
        seq_row = sequences_df[sequences_df['Protein Name'] == histone]
        if not seq_row.empty:
            sequences.append(seq_row['Sequence'].values[0])
    return sequences

def add_sequences_column(strategies_df, sequences_df):
    sequences_column = strategies_df['Histones'].apply(
        lambda histone_list: get_sequence_for_histones(histone_list, sequences_df)
    )
    histones_index = strategies_df.columns.get_loc('Histones')
    strategies_df.insert(histones_index + 1, 'Sequences', sequences_column)
    return strategies_df

def output_strategies_df(strategies_df):
    output_path = '../outputs/proteins_clustered_taxonomy_strategy_sequence_db.csv'
    strategies_df.to_csv(output_path, index=False)
    print(f"Updated strategies_df saved to {output_path}")

if __name__ == "__main__":
    main()