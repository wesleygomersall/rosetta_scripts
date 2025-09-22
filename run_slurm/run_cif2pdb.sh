#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=cif2pdb
#SBATCH --output=logs/cif2pdb-%x_%j.out
#SBATCH --error=logs/cif2pdb-%x_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00           
#SBATCH --mem=4G

# Create list of files for .cif -> .pdf conversion: `cif2pdb.txt`
FILE_LIST="cif2pdb.txt"

if ! [ -s "$FILE_LIST" ]; then
    echo "$FILE_LIST is empty. Exiting."
    exit
fi

SCRIPT="$HOME/rosetta_scripts/cif2pdb.py"

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate pymol

for file in $(cat $FILE_LIST); do
    python3 $SCRIPT -i $file
done

exit
