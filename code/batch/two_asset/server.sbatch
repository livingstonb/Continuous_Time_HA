#!/bin/bash
#SBATCH --job-name=ctime
#SBATCH --output=/home/livingstonb/GitHub/Continuous_Time_HA/output/two_asset/run%a.out
#SBATCH --error=/home/livingstonb/GitHub/Continuous_Time_HA/output/two_asset/run%a.err
#SBATCH --partition=broadwl
#SBATCH --array=14-31
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/Continuous_Time_HA/master_two_asset.m
