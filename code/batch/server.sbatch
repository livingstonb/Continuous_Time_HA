#!/bin/bash
#SBATCH --job-name=ctime
#SBATCH --output=/home/livingstonb/GitHub/Continuous_Time_HA/output/run%a.out
#SBATCH --error=/home/livingstonb/GitHub/Continuous_Time_HA/output/run%a.err
#SBATCH --partition=broadwl
#SBATCH --array=1,301
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/Continuous_Time_HA/master.m
