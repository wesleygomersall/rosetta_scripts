#!/usr/bin/env python3

# Wesley Gomersall
# Thu, 31 Jul 2025 17:51:27 -0700

import argparse
from pyrosetta import *

def cif_to_pdb(cifpath: str): 
    '''
    Creates pdb file from pose read from cif file.
    Input:
        cifpath: str        Path to cif file.
    '''
    pose = import_pose(cifpath)
    outputname = cifpath.strip(".cif")[0] + ".pdb"
    pose.dump_pdb(output)
    return

def main():
    parser = argparse.ArgumentParser(description="Create pdb file from cif.")
    parser.add_argument("input", "i", type=str, help=".cif file path")
    args = parser.parse_args()

    pyrosetta.init("-mute all")
    cif_to_pdb(args.input)

if __name__=="__main__":
    main()
