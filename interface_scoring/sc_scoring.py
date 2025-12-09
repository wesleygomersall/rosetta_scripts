#!/usr/bin/env python3 

import argparse
import numpy as np
import pyrosetta as pr
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector

parser = argparse.ArgumentParser(description="")
parser.add_argument("--file-list", "-f", type=str, help="Input pdb file list for interface scoring. One pdb path per line :-).")
parser.add_argument("--debug", action="store_true", default=False, help="Debug will output poses.") 
args = parser.parse_args()

pr.init()

def create_index_selector(res_nums):
    idx_selector = ResidueIndexSelector()
    for res in res_nums:
        idx_selector.append_index(res)
    return idx_selector

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
        for line in fin:
            file = line.strip()
            if file != '': my_file_list.append(file)

    for pdb_file in my_file_list: 
        with open(f"{pdb_file.rstrip('.pdb')}_residuescscore.csv", 'w') as fout:
            fout.write("Residue,sc_score\n")
            pose = pr.pose_from_pdb(pdb_file)
            
            binder_resis_beg = pose.conformation().chain_begin(2)
            binder_resis_end = pose.conformation().chain_end(2)

            for res in range(binder_resis_beg, binder_resis_end + 1):
                print(res)

                other_residues = list(range(1, pose.total_residue() + 1))
                other_residues.remove(int(res))

                if res < pose.total_residue() - 1:
                    sel2_string = f"{min(other_residues)}-{res-1},{res+1}-{pose.total_residue()}"
                elif res == pose.total_residue() - 1: 
                    sel2_string = f"{min(other_residues)}-{res-1},{pose.total_residue()}"
                elif res == pose.total_residue(): 
                    sel2_string = f"{min(other_residues)}-{res-1}"
                else: 
                    sel2_string = f"{min(other_residues)}-{pose.total_residue()}"

                sc = ShapeComplementarityFilter()
                sc.use_rosetta_radii(True)
                sc.residues1(str(res))
                sc.residues2(sel2_string) 
                sc_val = sc.score(pose)

                val = sc_val
                # val = myscores_all.sc_value
                fout.write(f"{res},{val}\n")

if __name__ == "__main__":
    main()
