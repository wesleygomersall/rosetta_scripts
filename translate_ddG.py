#!/usr/bin/env python3 

import argparse
import os
import datetime 
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
from pyrosetta.rosetta.protocols.relax import ClassicRelax
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
init("-mute all")

def unbind(pose, partners):
    """Perform translation movement to form unbound system

    PARAMS
    ------
        pose: pyrosetta.object
            bound complex pose
        partners: str
            Specifying the jump between A_B

    From Andrew Powers' 
    https://github.com/PowersPope/UsefulRosettaScripts/blob/main/scripts/calculate_ddg.py
    """
    # Set up the foldtree to allow unbinding
    setup_foldtree(pose, partners, Vector1([-1, -1, -1]))
    # Set up the translation mover and select what to move
    trans_mover = RigidBodyTransMover(pose, 1)
    # set number of steps to move
    trans_mover.step_size(100)
    trans_mover.apply(pose)

def main(): 
    parser = argparse.ArgumentParser(description="Calculate scores for relaxed protein-peptide complex and for translated peptide")
    parser.add_argument("--input", "-i", type=str, help="Path to input pdb file.")
    parser.add_argument("--refine", action="store_true", default=False, help="Relax structure, default will not.") 
    parser.add_argument("--dump", action="store_true", default=False, help="Dump pre and post-translated poses to pdb") 
    parser.add_argument("--design", action="store_true", default=False, help="Design new peptide sequences.") 
    parser.add_argument("--ndesigns", type=int, default=0, help="Number of designs.")
    args = parser.parse_args()
    
    if args.ndesigns > 0 and not args.design:
        print("Incompatible arguments: \
              Design sequences option not specified, \
              but number of designs not also zero. \
              Reassigning --ndesigns to zero.")
        args.ndesigns = 0 
    if args.ndesigns < 1 and args.design:
        print("Incompatible arguments: \
              Design sequences option specified, \
              but number of designs is less than 1. \
              Exiting.")
        return

    date_time_string = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    output_directory = date_time_string + "_ddGout_" + os.path.basename(args.input).split(".")[0] 
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    csv_out = open(output_directory + '/out_energies.csv', 'w')
    csv_out.write(f"DesignNum,PeptideSequence,bound_dG,unbound_dG,ddG\n")
    
    scorefxn = create_score_function("ref2015_cart")
    relax = ClassicRelax() 
    relax.set_scorefxn(scorefxn) 

    print(f"Input pose: {args.input}")
    input_pose = pose_from_pdb(args.input)
    input_score = scorefxn(input_pose)
    print(f"Input structure score: {input_score}")

    relaxed_input = input_pose.clone()
    if args.refine: 
        relax.apply(relaxed_input) 
        relaxed_score = scorefxn(relaxed_input)
        print(f"Relaxed score: {relaxed_score}")
    else: 
        print("Relax=False")
        relaxed_score = input_score

    for structnum in range(0, args.ndesigns + 1):
        print(f"Structure: {structnum}")
        
        test_pose = relaxed_input.clone()

        if structnum != 0 and args.design:
            # https://github.com/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.02-Packing-design-and-regional-relax.ipynb
            # https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.04-Protein-Design-2.ipynb
            chain_A = residue_selector.ChainSelector("A")   
            tf = TaskFactory()
            tf.push_back(operation.InitializeFromCommandline())
            tf.push_back(operation.IncludeCurrent())
            # Disable packing and design of chain A
            tf.push_back(operation.OperateOnResidueSubset( 
                operation.PreventRepackingRLT(), chain_A)) 

            # Convert to PackerTask to view
            # packer_task = tf.create_task_and_apply_taskoperations(test_pose)
            # if structnum == 1: print(packer_task)

            mm = MoveMap()
            mm.set_bb(True)
            mm.set_chi(True)
            mm.set_jump(True)

            rel_design = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1, script_file="MonomerDesign2019")
            rel_design.cartesian(True)
            rel_design.set_task_factory(tf)
            rel_design.set_movemap(mm)
            rel_design.minimize_bond_angles(True)
            rel_design.minimize_bond_lengths(True)
            rel_design.apply(test_pose) 
        
        peptide_sequence = test_pose.chain_sequence(2)
        print(f"Peptide: {peptide_sequence}") 
        # print(test_pose.pdb_info()) 

        bound_dG = scorefxn(test_pose)
        print(f"Bound score: {bound_dG}")

        if args.dump:
            pre_translated_out_path = f"{output_directory}/structure{structnum}_relaxed.pdb"
            test_pose.dump_pdb(pre_translated_out_path)
        
        unbind(test_pose, "A_B")

        if args.refine: 
            print("Relax transformed pose")
            relax.apply(test_pose) 

        unbound_dG = scorefxn(test_pose)
        print(f"Unbound score: {unbound_dG}")

        if args.dump: 
            translated_out_path = f"{output_directory}/structure{structnum}_translated.pdb"
            test_pose.dump_pdb(translated_out_path)

        ddG = bound_dG - unbound_dG
        print(f"ddG: {ddG}")

        csv_out.write(f"{structnum},{peptide_sequence},{bound_dG},{unbound_dG},{ddG}\n")

    csv_out.close()

if __name__ == "__main__": 
    main()
