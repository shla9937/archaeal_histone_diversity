#!/usr/bin/env python3

import pandas as pd
import argparse
import ast

def main():
    parser = argparse.ArgumentParser(description="Search for UniProt IDs based on histone sequences.")
    parser.add_argument('--df', type=str, required=True, help="Path to the CSV file with protein sequences.")
    args = parser.parse_args()

    sequences_df = pd.read_csv(args.df)
    uniprots = find_uniprot_ids(sequences_df)
    full_df = add_sequences_column(sequences_df, uniprots)
    output_strategies_df(full_df)

def find_uniprot_ids(sequences_df):
    # sequences = []
    # for histone in histone_list:
    #     seq_row = sequences_df[sequences_df['Protein Name'] == histone]
    #     if not seq_row.empty:
    #         sequences.append(seq_row['Sequence'].values[0])
    # return sequences
    print(sequences_df['Sequences'])
    True

def add_sequences_column(strategies_df, sequences_df):
    # sequences_column = strategies_df['Histones'].apply(
    #     lambda histone_list: get_sequence_for_histones(histone_list, sequences_df)
    # )
    # histones_index = strategies_df.columns.get_loc('Histones')
    # strategies_df.insert(histones_index + 1, 'Sequences', sequences_column)
    # return strategies_df
    True

def output_strategies_df(strategies_df):
    # output_path = '../outputs/proteins_clustered_taxonomy_strategy_sequence_db.csv'
    # strategies_df.to_csv(output_path, index=False)
    # print(f"Updated strategies_df saved to {output_path}")
    True
if __name__ == "__main__":
    main()