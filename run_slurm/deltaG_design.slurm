#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=rd_dG
#SBATCH --output=logs/pyrosetta_dG_design-%x_%j.out
#SBATCH --error=logs/pyrosetta_dG_design-%x_%j.err
#SBATCH --partition=computelong
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00           
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --mail-user=wesg@uoregon.edu

if [ $# != 2 ]; then
    echo "Supply 2 arguments: path/to/input.pdb numdesigns"
    exit
fi

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate pyrosetta

python3 $HOME/rosetta_scripts/translate_ddG.py \
    --input $1 \
    --refine \
    --dump \
    --design \
    --nstructs $2
