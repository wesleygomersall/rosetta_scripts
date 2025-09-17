#!/usr/bin/env python3 

import argparse
from pyrosetta import init
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.protocols.relax import ClassicRelax

def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, help="Path to input pdb file.")
    args = parser.parse_args()

    scorefxn = 'ref2015'
    relaxed_out_path = "relaxed.pdb"
    # translated_out_path = "translated.pdb" 

    init("-mute all")

    pose = pose_from_pdb(args.input)
    
    relax = ClassicRelax() 
    relax.set_scorefxn(scorefxn) 
    relax.apply(pose) 

    # TODO 
    # score delta G 

    pose.dump_pdb(relaxed_out_path)

    # TODO 
    # translate peptide from translated pose
        # first find vector between two chains' centers of mass
        # translate a chain relative to another chain  from PyRosetta manual:
            # trans_mover = RigidBodyTransMover( pose, jump_num )
            # trans_mover.step_size(50)
            # trans_mover.apply( pose )
        # relax -maybe
        # Export pdb from translated pose 
        # score delta G 

    # TODO 
    # calc delta delta G 

if __name__ == "__main__": 
    main()
