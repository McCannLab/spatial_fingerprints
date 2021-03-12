#!/bin/bash
#SBATCH --account=def-ksmccann
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-10
#SBATCH --mail-user=kcazelle@uoguelph.ca
#SBATCH --mail-type=FAIL

julia -p 1 main_all123.jl 1 $SLURM_ARRAY_TASK_ID 20 10 10 700 0