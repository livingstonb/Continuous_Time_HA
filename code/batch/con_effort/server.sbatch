#!/bin/bash
#SBATCH --job-name=coneffort
#SBATCH --output=/home/livingstonb/GitHub/Continuous_ConEffort/output/run%a.out
#SBATCH --error=/home/livingstonb/GitHub/Continuous_ConEffort/output/run%a.err
#SBATCH --partition=broadwl
#SBATCH --array=1-32
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/Continuous_ConEffort/master.m
