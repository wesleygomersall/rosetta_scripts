#!/usr/bin/env python3

# Wesley Gomersall 
# May 27th 2025 
# University of Oregon, P.I.: Parisa Hosseinzadeh

import os

from pyrosetta import * 
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory
from pyrosetta.rosetta.protocols.docking import setup_foldtree

# options 
# input_pdb= "/home/wesg/rosetta_fpd/comd/comd-c16_extended.pdb" # use absolute path
# The above is the extended conformation peptide structure
input_pdb = "/home/wesg/pyrosetta_scripts/input.pdb"
output = "/home/wesg/pyrosetta_scripts/comdn230-c16_output_testing"
partners = "A_B"
ligand_params = "" # need params for non-peptide small molecule docking
jobs = 10 # number of trajectories
job_output = "lig_dock_output" 

# dir to store the outputs 
if not os.path.isdir(output):
    os.mkdir(output)
os.chdir(output) 

def main():
    """
    Perform ligand-protein docking with fullatom docking (DockingHighRes) on
    complex in input PDB. 
    Adapted from Rosetta/main/source/src/python/PyRosetta/src/demo/D120_Ligand_interface.py
    
    Uses default params, and works only for two body docking between peptide 
    chains! 
    """
    pyrosetta.init("-mute all") # required

    # 1. set up input pose
    pose = pose_from_pdb(input_pdb)
    original_pose = pose.clone()
    
    # 2. set up pose FoldTree for docking 
    # Just using default FoldTree with single jump between chains
    dock_jump = 1
    # pyrosetta.rosetta.protocols.docking.setup_foldtree(pose, partners, Vector1([dock_jump]))
    # print(pose.fold_tree())
    # pyrosetta.rosetta.protocols.docking.set_autofoldtree(True)
    setup_foldtree(pose, partners, Vector1([dock_jump]))
    
    # 3. create copy of pose to be modified
    test_pose = Pose()
    test_pose.assign(pose)

    # 4. create score functions for scoring ligand-protein complexes
    sfxn = create_score_function('ligand')

    # 5. set up the DockMCMProtocol object for fullatom docking
    docking = pyrosetta.rosetta.protocols.docking.DockMCMProtocol()
    docking.set_scorefxn(sfxn)

    # 6. create PyJobDistributor for managing trajectories
    jd = PyJobDistributor(job_output, jobs, sfxn)

    # 7. (optional) create PyMOL_Observer for viewing intermediate output

    # 9. perform protein-protein docking
    counter = 0 
    while not jd.job_complete: 
        test_pose.assign(pose)
        counter += 1
        test_pose.pdb_info().name(job_output + '_' + str(counter))

        docking.apply(test_pose)

        jd.output_decoy(test_pose)

if __name__=="__main__":
    # testing()
    main()
