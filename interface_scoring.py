#!/usr/bin/env python3 

import argparse
# import os
# import datetime 
# import pandas as pd
import pyrosetta as pr
# from pyrosetta.rosetta.core.scoring import *
# from pyrosetta.rosetta.core.select import residue_selector
# from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", "-i", type=str, help="Input pdb file for interface scoring.")
parser.add_argument("--debug", action="store_true", default=False, help="Debug will output poses.") 
args = parser.parse_args()

print(args.input)

pr.init()

def sc_score_per_res(input_pose, 
                     target_chainid = 1, 
                     binder_chainid = 2,
                     debug = False):

    target_resis = input_pose.chain_sequence(target_chainid)
    binder_resis = input_pose.chain_sequence(binder_chainid)
    interface_scores = []

    for res_num in range(1, len(binder_resis) + 1):
        test_pose = pr.Pose()
        test_pose.assign(input_pose)

        # delete all resis not equal to current residue number from binder chain
        for i in range(len(binder_resis)):
            if len(target_resis) + len(binder_resis) - i != res_num + len(target_resis):
                test_pose.delete_residue_slow(len(target_resis) + len(binder_resis) - i)

        # score interface between single residue and target
        interface = f"A_B:{res_num}" # not for input to score_interface
        ifs = score_interface(test_pose, "A_B")
        interface_scores.append([res_num, binder_resis[res_num - 1], interface, ifs.sc_value])

        if debug: test_pose.dump_pdb(f"test_{interface}.pdb")

    return interface_scores

def score_interface(pose, interface ):
    # analyze interface statistics
    iam = InterfaceAnalyzerMover()
    iam.set_interface(interface)
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_compute_interface_delta_hbond_unsat(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
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
    pdb_file = args.input
    pose = pr.pose_from_pdb(pdb_file)

    myscores = sc_score_per_res(pose, debug = args.debug)

    with open(f"{args.input.strip('.pdb')}_sc_per_res.csv", 'w') as fout:
        fout.write(f"Residue_number,Residue,Interface,sc_value\n")
        for item in myscores:
            fout.write(','.join(map(str, item)) + '\n')

if __name__ == "__main__":
    main()
