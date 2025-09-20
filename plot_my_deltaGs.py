#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

file_list = ["ComD_C16_AF3_allalanine_fixed_dgoutput/out_energies.csv",
             "ComD_C16_AF3_fixed_dgoutput/out_energies.csv",
             "comd_c16g2_AF3_fixed_dgoutput/out_energies.csv",
             "comd_m8g2_AF3_fixed_dgoutput/out_energies.csv",
             "eckertC16-11_fixed_dgoutput/out_energies.csv",
             "eckertC16-12_fixed_dgoutput/out_energies.csv",
             "eckertC16-13_fixed_dgoutput/out_energies.csv",
             "eckertC16-14_fixed_dgoutput/out_energies.csv",
             "eckertC16-15_fixed_dgoutput/out_energies.csv",
             "eckertC16-16_fixed_dgoutput/out_energies.csv",
             "eckertC16-17_fixed_dgoutput/out_energies.csv",
             "eckertC16-18_fixed_dgoutput/out_energies.csv",
             "mpnn_neg_ctrl1_fixed_dgoutput/out_energies.csv",
             "sepmcsp21_model_fixed_dgoutput/out_energies.csv"]

def collect_to_one_file(file_list: list):
    with open("deltaGs.csv", 'w') as fout: 
        fout.write("structure,deltaG\n")
        for file in file_list: 
            df = pd.read_csv(file)
            name = file.split("_fixed_dgoutput")[0]
            fout.write(f"{name},{df.iloc[3, 1]}\n")

def plot_dGs(filename: str): 
    plt.clf()
    df = pd.read_csv(filename)
    fig, ax = plt.subplots()
    labels = df['structure']
    ax.bar(range(len(labels)),df['deltaG'])

    ax.set_xlabel('Variant')
    ax.set_xticks(range(len(labels)), labels, rotation='vertical')
    ax.set_ylabel('dG')
    ax.set_title('Delta G: CSP variants')
    plt.subplots_adjust(bottom=0.4) 
    plt.savefig("dG_figure.pdf", format="pdf")

if __name__ == "__main__":
    plot_dGs("deltaGs.csv")
    
    

