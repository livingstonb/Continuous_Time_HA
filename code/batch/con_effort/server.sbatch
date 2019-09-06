#!/bin/bash
#SBATCH --job-name=coneffort
#SBATCH --output=/home/livingstonb/GitHub/Continuous_Time_HA/output/con_effort/run%a.out
#SBATCH --error=/home/livingstonb/GitHub/Continuous_Time_HA/output/con_effort/run%a.err
#SBATCH --partition=broadwl
#SBATCH --array=1-36
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load matlab

matlab -nodisplay < /home/livingstonb/GitHub/Continuous_Time_HA/master_con_effort.m
