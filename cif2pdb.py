#!/usr/bin/env python3

# Wesley Gomersall
# Mon, 22 Sep 2025 13:52:49 -0700

import argparse
import pymol2

def cif_to_pdb(cifpath: str): 
    '''
    Creates pdb file from pose read from cif file. Requires PyMOL.
    Input:
        cifpath: str        Path to cif file.
    '''
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(cifpath, 'mystructure')
        pymol.cmd.save(cifpath.replace('.cif', '.pdb'), selection='mystructure')
    return

def main():
    parser = argparse.ArgumentParser(description="Create pdb file from cif.")
    parser.add_argument("--input", "-i", type=str, help=".cif file path")
    args = parser.parse_args()
    cif_to_pdb(args.input)

if __name__=="__main__":
    main()
