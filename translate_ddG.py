#!/usr/bin/env python3 

import argparse
import os
import datetime 
import pandas as pd
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
from pyrosetta.rosetta.protocols.relax import ClassicRelax
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.protocols.residue_selectors import StoredResidueSubsetSelector, StoreResidueSubsetMover

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
    parser.add_argument("--refine", action="store_true", default=False, help="Relax input structure, default will not.") 
    parser.add_argument("--dump", action="store_true", default=False, help="Dump pre and post-translated poses to pdb") 
    parser.add_argument("--design", action="store_true", default=False, help="Design new peptide sequences.") 
    parser.add_argument("--nstructs", type=int, default=0, help="Number of designs to perform. If not using --design option, this adds replicates to the translate step.")
    parser.add_argument("--debug", action="store_true", default=False, help="Debug mode, dont apply relax/design instead print packers") 
    args = parser.parse_args()
    
    if args.nstructs < 1 and args.design:
        print("Incompatible arguments: \
              Design sequences option specified, \
              but number of designs is less than 1. \
              Exiting.")
        return

    date_time_string = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    date_string = datetime.datetime.now().strftime("%Y%m%d")
    output_parent_dirname = "output/" + date_string + "_deltaGoutput/"
    output_sub_dirname = date_time_string + "_ddGout_" + os.path.basename(args.input).split(".")[0] 

    output_directory = output_parent_dirname + output_sub_dirname
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    scorefxn = create_score_function("ref2015_cart")
    relax = ClassicRelax() 
    relax.set_scorefxn(scorefxn) 

    print(f"Input pose: {args.input}")
    input_pose = pose_from_pdb(args.input)
    input_score = scorefxn(input_pose)
    print(f"Input structure score: {input_score}")

    relaxed_input = input_pose.clone()
    if args.refine and not args.debug: 
        relax.apply(relaxed_input) 
        relaxed_score = scorefxn(relaxed_input)
        print(f"Relaxed score: {relaxed_score}")
    else: 
        print("Relax=False")
        relaxed_score = input_score

    all_structnum = []
    all_peptide_sequence = []
    all_bound_dG = []
    all_unbound_dG = []
    all_ddG = []

    # https://github.com/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.02-Packing-design-and-regional-relax.ipynb
    # https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.04-Protein-Design-2.ipynb
    tf = TaskFactory()
    tf.clear()
    tf.push_back(operation.InitializeFromCommandline())

    chA_selector = residue_selector.ChainSelector("A")   
    chB_selector = residue_selector.NotResidueSelector(chA_selector)

    nbr_selector = residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(chB_selector)
    nbr_selector.set_include_focus_in_subset(False)

    interface_selector = residue_selector.OrResidueSelector()
    interface_selector.add_residue_selector(nbr_selector)
    interface_selector.add_residue_selector(chB_selector)

    stored_interface = StoredResidueSubsetSelector("neighbor_pre_translation")

    # Generate selection storage & apply
    store_mover = StoreResidueSubsetMover(interface_selector,
                                          "neighbor_pre_translation",
                                          True
                                          )
    store_mover.apply(relaxed_input)

    # prevent repacking on everything outside of nbr_selector & chB
    tf.push_back(operation.OperateOnResidueSubset(
        operation.PreventRepackingRLT(), stored_interface, True ))

    # for chain A neighbors of chain B, repack only
    tf.push_back(operation.OperateOnResidueSubset( 
        operation.RestrictToRepackingRLT(), nbr_selector)) 

    for structnum in range(0, args.nstructs + 1):
        print(f"\nStructure: {structnum}")
        all_structnum.append(structnum)
        
        test_pose = relaxed_input.clone()
        
        tf_new = tf.clone()

        if structnum == 0 or not args.design: # also restrict Chain B to repack only
            tf_new.push_back(operation.OperateOnResidueSubset( 
                operation.RestrictToRepackingRLT(), chB_selector)) 
            
        if args.debug: 
            # Convert tf to PackerTask to view
            print("before move")
            packer_task = tf_new.create_task_and_apply_taskoperations(test_pose)
            print(packer_task)

        mm = MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)
        mm.set_jump(True)

        rel_design = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1, script_file="MonomerDesign2019")
        rel_design.cartesian(True)
        rel_design.set_task_factory(tf_new)
        rel_design.set_movemap(mm)
        rel_design.minimize_bond_angles(True)
        rel_design.minimize_bond_lengths(True)
        if not args.debug: 
            rel_design.apply(test_pose) 
        
        all_peptide_sequence.append(test_pose.chain_sequence(2))
        print(f"Peptide: {all_peptide_sequence[-1]}") 
        # print(test_pose.pdb_info()) 

        all_bound_dG.append(scorefxn(test_pose))
        print(f"Bound score: {all_bound_dG[-1]}")

        if args.dump:
            pre_translated_out_path = f"{output_directory}/structure{structnum}_relaxed.pdb"
            test_pose.dump_pdb(pre_translated_out_path)
        
        unbind(test_pose, "A_B")

        # Always want to avoid re-designing sequence after unbinding
        tf_new.push_back(operation.OperateOnResidueSubset( 
            operation.RestrictToRepackingRLT(), stored_interface)) 
        rel_design.set_task_factory(tf_new)

        if args.debug: 
            print("after move")
            packer_task = tf_new.create_task_and_apply_taskoperations(test_pose)
            print(packer_task)

        if not args.debug: 
            # Re-relax without design
            rel_design.apply(test_pose) 

        all_unbound_dG.append(scorefxn(test_pose))
        print(f"Unbound score: {all_unbound_dG[-1]}")

        if args.dump: 
            translated_out_path = f"{output_directory}/structure{structnum}_translated.pdb"
            test_pose.dump_pdb(translated_out_path)

        all_ddG.append(all_bound_dG[-1] - all_unbound_dG[-1])
        print(f"ddG: {all_ddG[-1]}")

    allstructs_ddG_df = pd.DataFrame({
        "DesignNum": all_structnum,
        "PeptideSequence": all_peptide_sequence,
        "bound_dG": all_bound_dG,
        "unbound_dG": all_unbound_dG,
        "ddG": all_ddG})
    allstructs_ddG_df.to_csv(f"{output_directory}/out_energies.csv", index=False)

if __name__ == "__main__": 
    main()
