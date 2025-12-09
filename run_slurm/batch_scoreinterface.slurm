#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=ifsc
#SBATCH --output=logs/interface_scoring-%x_%j.out
#SBATCH --error=logs/interface-scoring%x_%j.err
#SBATCH --array=0-10                # Job array for 10 tasks (adjust as needed)
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00           
#SBATCH --mem=8G
#SBATCH --mail-type=END
#SBATCH --mail-user=wesg@uoregon.edu 

SC_SCRIPT="$HOME/rosetta_scripts/interface_scoring/sc_scoring.py"
PS_SCRIPT="$HOME/rosetta_scripts/interface_scoring/packstat_perres.py"

# Create list of files for .cif -> .pdf conversion: `inputs_cif2pdb.txt`
FILE_LIST="inputs_scoreinterface.txt"

if ! [ -s "$FILE_LIST" ]; then
    echo "$FILE_LIST is empty. Exiting."
    exit
fi

total_files=$(wc -l < $FILE_LIST)
files_per_task=$(( (total_files + SLURM_ARRAY_TASK_COUNT - 1) / SLURM_ARRAY_TASK_COUNT ))  # Round up
start_index=$(( SLURM_ARRAY_TASK_ID * files_per_task ))
end_index=$(( start_index + files_per_task - 1 ))

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate pyrosetta

echo "Task ID: $SLURM_ARRAY_TASK_ID, Processing files $start_index to $end_index"

sed -n "$((start_index + 1)),$((end_index + 1))p" $FILE_LIST | while read file_path; do
    python3 $PS_SCRIPT -i ${file_path}
    python3 $SC_SCRIPT -i ${file_path}
done

exit
