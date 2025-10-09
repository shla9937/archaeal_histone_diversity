#!/usr/bin/env python

import pandas as pd
import argparse
import ast

def main():
    parser = argparse.ArgumentParser(description="Search for UniProt IDs based on histone sequences.")
    parser.add_argument('--df', type=str, required=True, help="Path to the CSV file with protein sequences.")
    parser.add_argument('--ids', type=str, required=True, help="Path to the TSV file with accessions.")
    args = parser.parse_args()

    protein_df = pd.read_csv(args.df)
    uniprots = pd.read_csv(args.ids,sep='\t',header=None,names=['UniProtID', 'Score'],dtype=str)
    full_df = add_sequences_column(protein_df, uniprots)
    output_strategies_df(full_df)

def add_sequences_column(protein_df, uniprots):
    protein_df.insert(4, 'UniProtIDs', None)
    protein_df['Histones'] = protein_df['Histones'].apply(ast.literal_eval)
    for index, row in protein_df.iterrows():
        uniprot_list = []
        for histone in row['Histones']:
            for index_id, row_id in uniprots.iterrows():
                if index_id == histone.lstrip('>'):
                    uniprot_list.append(row_id.UniProtID)
        protein_df.loc[index, 'UniProtIDs'] = str(uniprot_list)
    return protein_df

def output_strategies_df(strategies_df):
    output_path = '../proteins_clustered_taxonomy_strategy_sequence_uniprot_db.csv'
    strategies_df.to_csv(output_path, index=False)

if __name__ == "__main__":
    main()