#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
from decimal import Decimal, getcontext

def main():
    parser = argparse.ArgumentParser(description='Plot data from multiple files.')
    parser.add_argument('--files', nargs='+', type=str, required=True, help='List of input files')
    parser.add_argument('--names', nargs='+', type=str, required=True, help='List of names corresponding to each file')
    args = parser.parse_args()

    if len(args.files) != len(args.names):
        raise ValueError("The number of files and names must be the same.")

    global cud_palette
    cud_palette = ["#56B4E9","#D55E00","#E69F00","#F0E442","#009E73","#CC79A7","#0072B2","#000000", "#999999"]

    getcontext().prec = 100

    plt.figure(figsize=(2.5, 2))
    plt.rcParams.update({'font.size': 6})

    for i, file in enumerate(args.files):
        data = pd.read_csv(file, sep="\t", header=None, dtype=str)
        data[0] = data[0].apply(Decimal)
        data[1] = data[1].apply(Decimal)
        plt.plot(data[0], data[1], label=args.names[i], color=cud_palette[i % len(cud_palette)])


    plt.xscale('log')
    plt.xlim(1e-40, 1e10)
    plt.xlabel("E-value")
    plt.ylabel("Hits")
    # plt.title("Hits at a give E-value")
    plt.legend(loc='best')
    plt.tight_layout()
    
    plt.savefig("e_value_comparison.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()