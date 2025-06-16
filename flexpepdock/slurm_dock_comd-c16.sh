#!/usr/bin/env bash

#SBATCH --job-name=pros
#SBATCH --output=/home/wesg/pyrosetta_scripts/slurm_logstd.out
#SBATCH --error=/home/wesg/pyrosetta_scripts/slurm_logstd.err

#SBATCH --account=parisahlab
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

#SBATCH --mail-user=wesg@uoregon.edu
#SBATCH --mail-type=END

conda activate pyrosetta

cd ~/pyrosetta_scripts

python3 dock_comd-c16.py

conda deactivate
