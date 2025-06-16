#!/usr/bin/env python3

# Wesley Gomersall
# May 28, 2025

# Input two pdb files and calculate the RMSD between the alpha carbons of two 
# poses.

import argparse

from pyrosetta import *
from pyrosetta.rosetta.core.scoring import CA_rmsd

def print_chains_seqs(mypose):
    """
    For input pose, print length and sequences of each chain.
    """
    # Get number of chains in structure
    numchains = mypose.chain(len(mypose.sequence()))
    for i in range(numchains):
        print(f"Chain {i+1}:")
        seq = mypose.chain_sequence(i+1)
        print(f"    Length: {len(seq)}")
        print(f"    Sequence: {seq}")

def main():
    parser = argparse.ArgumentParser(description="Calculate the RMSD between alpha carbons of two structures.")
    parser.add_argument("pose1", type=str, help="The first structure.")
    parser.add_argument("pose2", type=str, help="The second structure.")
    args = parser.parse_args()

    pyrosetta.init("-mute all")

    first_pose = pose_from_pdb(args.pose1)
    seq_len1 = len(first_pose.sequence())
    second_pose = pose_from_pdb(args.pose2)
    seq_len2 = len(second_pose.sequence())
    ca_rmsd = CA_rmsd(first_pose, second_pose)
    
    print("whole-structure RMSD:")
    print(f"\nStructure 1: {args.pose1}")
    print_chains_seqs(first_pose)
    print(f"\nStructure 2: {args.pose2}")
    print_chains_seqs(second_pose)
    print(f"\nRMSD = {ca_rmsd} ({seq_len1} to {seq_len2} alpha carbons)")

if __name__=="__main__":
    main()
