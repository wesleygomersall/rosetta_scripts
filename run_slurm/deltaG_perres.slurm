#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=pres_ddG
#SBATCH --output=logs/pyrosetta_dG_design-%x_%j.out
#SBATCH --error=logs/pyrosetta_dG_design-%x_%j.err
#SBATCH --partition=computelong
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00           
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --mail-user=wesg@uoregon.edu

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate pyrosetta

FILE_LIST="inputs_perres_ddG.txt"
if ! [ -s "$FILE_LIST" ]; then
    echo "$FILE_LIST is empty. Exiting."
    exit
fi

for item in $(cat $FILE_LIST); do

    python3 ~/rosetta_scripts/translate_ddG.py \
        -i "$item" \
        --refine \
        --dump \
        --perres-ddg \
        # --debug

done

conda deactivate
