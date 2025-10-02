#!/bin/bash

# This script initiates jobs with `deltaG_complex.slurm` for each input pdb file
# in `inputs_dGcomplex.txt`

FILE_LIST="inputs_dGcomplex.txt"
if ! [ -s "$FILE_LIST" ]; then
    echo "$FILE_LIST is empty. Exiting."
    exit
fi

SLURM_SCRIPT="$HOME/rosetta_scripts/run_slurm/deltaG_complex.slurm"

for file in $(cat $FILE_LIST); do
    sbatch $SLURM_SCRIPT $file
done

exit
