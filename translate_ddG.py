#!/usr/bin/env python3 

import argparse
import os
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.relax import ClassicRelax
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree

def main(): 
    parser = argparse.ArgumentParser(description="Calculate scores for relaxed protein-peptide complex and for translated peptide")
    parser.add_argument("--input", "-i", type=str, help="Path to input pdb file.")
    parser.add_argument("--refine", action="store_true", default=False, help="Relax structure, default will not.") 
    args = parser.parse_args()

    init("-mute all")

    scorefxn = create_score_function("ref2015") 
    # print(scorefxn) 
    relax = ClassicRelax() 
    relax.set_scorefxn(scorefxn) 

    output_directory = os.path.basename(args.input).split(".")[0] + "_dgoutput"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    pre_translated_out_path = output_directory + "/structure_relaxed.pdb"
    translated_out_path = output_directory + "/structure_translated.pdb" 

    print(f"Input pose {args.input}")
    input_pose = pose_from_pdb(args.input)
    input_score = scorefxn(input_pose)
    print(f"input structure score: {input_score}")
    relaxed_pose = Pose()
    relaxed_pose.assign(input_pose)

    if args.refine: # relax and re score
        relax.apply(relaxed_pose) 
        relaxed_score = scorefxn(relaxed_pose)
        print(f"relaxed score: {relaxed_score}")
    else: 
        print("Relax=False")
        relaxed_score = input_score

    relaxed_pose.dump_pdb(pre_translated_out_path)
    translated_pose = Pose()
    translated_pose.assign(relaxed_pose)

    # Translate chain B (peptide) away from chain A (protein)
    setup_foldtree(translated_pose, "A_B", Vector1([-1, -1, -1]))
    trans_mover = RigidBodyTransMover(translated_pose, 1) 
    trans_mover.step_size(50)
    trans_mover.apply(translated_pose)

    if args.refine: 
        print("relax transformed pose")
        relax.apply(translated_pose) 

    translated_score = scorefxn(translated_pose)
    print(f"post-translation score: {translated_score}")
    translated_pose.dump_pdb(translated_out_path)

    # dG_binding = dG_bound - (dG_unboundpeptide + dG_unboundprotein)
    delta_g = relaxed_score - translated_score
    print(delta_g)

    with open(output_directory + '/out_energies.csv', 'w') as fout: 
        fout.write("struct,energy\n")
        fout.write(f"{args.input}_RAW,{input_score}\n")
        fout.write(f"{args.input}_relaxed,{relaxed_score}\n")
        fout.write(f"{args.input}_translated,{translated_score}\n")
        fout.write(f"{args.input}_diff,{delta_g}\n")

if __name__ == "__main__": 
    main()
