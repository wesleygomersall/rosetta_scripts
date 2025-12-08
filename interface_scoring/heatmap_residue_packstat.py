#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse 

def main():

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", nargs='+', type=str, default = "", help="")
    args = parser.parse_args()

    packstats = pd.DataFrame()
    for file in args.input: 
        name = file.split('/')[-1].rstrip('_residuepackstat.csv')
        temp = pd.read_csv(file)
        packstats[name] = temp.loc[:, 'packstat_value']

    plt.figure(figsize=(10, 5)) # Adjust figure size
    sns.heatmap(packstats) 
    plt.title("packstat per residue")
    plt.ylabel("Peptide residue number")
    plt.subplots_adjust(bottom=0.4) 
    plt.show()

if __name__ == "__main__": 
    main()

