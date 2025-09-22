#!/usr/bin/env python3

# Wesley Gomersall
# Thu, 31 Jul 2025 17:51:27 -0700

import argparse
import pyrosetta
from pyrosetta import dump_pdb
import pyrosetta.rosetta.core.import_pose as import_pose

def cif_to_pdb(cifpath: str): 
    '''
    Creates pdb file from pose read from cif file.
    Input:
        cifpath: str        Path to cif file.
    '''
    pose = import_pose.pose_from_file(cifpath)
    outputname = cifpath.strip(".cif") + ".pdb"
    print(f"writing pdb to {outputname}")
    pose.dump_pdb(outputname)
    return

def main():
    parser = argparse.ArgumentParser(description="Create pdb file from cif.")
    parser.add_argument("--input", "-i", type=str, help=".cif file path")
    args = parser.parse_args()

    pyrosetta.init("-mute all")
    cif_to_pdb(args.input)

if __name__=="__main__":
    main()
