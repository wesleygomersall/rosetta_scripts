#!/usr/bin/env python3 

import argparse
import numpy as np
import pyrosetta as pr
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

parser = argparse.ArgumentParser(description="")
parser.add_argument("--file-list", "-f", type=str, help="Input pdb file list for interface scoring. One pdb path per line :-).")
parser.add_argument("--debug", action="store_true", default=False, help="Debug will output poses.") 
args = parser.parse_args()

pr.init()

def score_interface(pose, interface ):
    # analyze interface statistics
    iam = InterfaceAnalyzerMover()
    iam.set_interface(interface)
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(False)
    iam.set_compute_interface_delta_hbond_unsat(False)
    iam.set_calc_dSASA(False)
    iam.set_calc_hbond_sasaE(False)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(False)
    iam.apply(pose)

    iam.apply(pose)
    
    # retrieve statistics
    interfacescore = iam.get_all_data()
    # interface_sc = interfacescore.sc_value # shape complementarity
    # interface_nres = iam.get_num_interface_residues() # number of interface residues
    # interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
    # interface_dG = iam.get_interface_dG() # interface dG
    # interface_dSASA = iam.get_interface_delta_sasa() # interface dSASA (interface surface area)
    # interface_packstat = iam.get_interface_packstat() # interface pack stat score

    return interfacescore

def main():

    # get list of inputs
    my_file_list = []
    with open(args.file_list, 'r') as fin: 
        line = fin.read().strip()
        if line != '': my_file_list.append(line)

    for pdb_file in my_file_list: 
        pose = pr.pose_from_pdb(pdb_file)

        myscores_all = score_interface(pose, "A_B")

        print(f"file: {pdb_file}")
        print(f"interface sc_value: {myscores_all.sc_value}")

if __name__ == "__main__":
    main()
