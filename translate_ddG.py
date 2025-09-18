#!/usr/bin/env python3 

import argparse
# import numpy as np
from pyrosetta import *

from pyrosetta.rosetta.core.scoring import *
# from pyrosetta.rosetta.core.pack.task import *
# from pyrosetta.rosetta.protocols import *
# from pyrosetta.rosetta.protocols.geometry import *

# from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.relax import ClassicRelax

def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, help="Path to input pdb file.")
    args = parser.parse_args()

    scorefxn = create_score_function("ref2015")
    print(scorefxn) 
    relaxed_out_path = "relaxed.pdb"
    translated_out_path = "translated.pdb" 

    init("-mute all")

    # Read pose get initial score
    in_pose = pose_from_pdb(args.input)
    pre_relaxed_score = scorefxn(in_pose)
    print(f"pre-relax score: {pre_relaxed_score}")
    
    # Relax and score
    relaxed_pose = Pose()
    relaxed_pose.assign(in_pose)
    relax = ClassicRelax() 
    relax.set_scorefxn(scorefxn) 
    relax.apply(relaxed_pose) 
    relaxed_score = scorefxn(relaxed_pose)
    print(f"post-relax score: {relaxed_score}")
    relaxed_pose.dump_pdb(relaxed_out_path)

    # Translate and score
    translated_pose = Pose()
    translated_pose.assign(relaxed_pose)

    # TODO: translate peptide from relaxed pose
    # com_chain0 = center_of_mass(select(chain0))
    # com_chain1 = center_of_mass(select(chain1))
    # direction = com_chain0 - com_chain1
    # direction = direction - np.linalg.norm(direction) 
    # translate a chain relative to another chain  from PyRosetta manual:
        # trans_mover = RigidBodyTransMover( pose, jump_num )
        # trans_mover.step_size(50)
        # trans_mover.apply( pose )

    translated_score = scorefxn(translated_pose)
    print(f"post-translation score: {translated_score}")
    translated_pose.dump_pdb(translated_out_path)

    # Get difference 
    delta_g = relaxed_score - translated_score
    print(delta_g)

if __name__ == "__main__": 
    main()
