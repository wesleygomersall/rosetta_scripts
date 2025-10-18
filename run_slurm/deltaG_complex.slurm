#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=deltaG_TransAB
#SBATCH --output=logs/pyrosetta_dG_trans-%x_%j.out
#SBATCH --error=logs/pyrosetta_dG_trans-%x_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00           
#SBATCH --mem=16G

if [ $# != 1 ]; then
    echo "Supply 2 arguments: path/to/input.pdb nstructs"
    exit
fi

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate pyrosetta

python3 $HOME/rosetta_scripts/translate_ddG.py \
    --input $1 \
    --refine \
    --dump \
    --nstructs $2
