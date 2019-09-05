#!/bin/bash
#SBATCH --job-name=ce_table
#SBATCH --output=/home/livingstonb/GitHub/Continuous_ConEffort/output/create_table.out
#SBATCH --error=/home/livingstonb/GitHub/Continuous_ConEffort/output/create_table.err
#SBATCH --partition=broadwl
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8000

module load matlab

matlab -nodisplay < "/home/livingstonb/GitHub/Continuous_ConEffort/code/batch/create_table.m"
