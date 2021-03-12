#!/bin/bash
#SBATCH --account=def-ksmccann
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-17
#SBATCH --mail-user=kcazelle@uoguelph.ca
#SBATCH --mail-type=FAIL

julia -p 1 main.jl 0 $SLURM_ARRAY_TASK_ID 20 10 10000 500 1