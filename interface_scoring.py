#!/usr/bin/env python3 

import argparse
import pyrosetta as pr
from pyrosetta.rosetta.core.scoring.packstat import pose_to_pack_data, compute_residue_packing_score

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", "-i", type=str, help="Input pdb file for interface scoring.")
parser.add_argument("--debug", action="store_true", default=False, help="Debug will output poses.") 
args = parser.parse_args()

pr.init()

def packstat_per_res(input_pose, 
                     target_chainid = 1,
                     binder_chainid = 2):

    target_resis = input_pose.chain_sequence(target_chainid)
    binder_resis = input_pose.chain_sequence(binder_chainid)
    start_index = len(target_resis) + 1
    end_index = start_index + len(binder_resis)

    packdata = pose_to_pack_data(input_pose)
    print()

    interface_scores = []

    for res_num in range(start_index, end_index):
        interface_scores.append((res_num, 
                                 input_pose.residue(res_num).name3(), 
                                 compute_residue_packing_score(packdata, res_num),
                                 ))


    return interface_scores

def main():
    pdb_file = args.input
    pose = pr.pose_from_pdb(pdb_file)

    myscores = packstat_per_res(pose)

    with open(f"{args.input.strip('.pdb')}_residuepackstat.csv", 'w') as fout:
        fout.write(f"Res_in_pdb,Res_in_chain,Residue,packstat_value\n")
        for i, item in enumerate(myscores):
            fout.write(f"{item[0]},{i+1},{item[1]},{item[2]}\n")

if __name__ == "__main__":
    main()
