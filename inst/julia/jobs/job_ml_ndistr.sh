#!/bin/bash
#SBATCH --account=def-ksmccann
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-200
#SBATCH --mail-user=kcazelle@uoguelph.ca
#SBATCH --mail-type=FAIL

nbios=(1 2 5 10 15)
nbio=${nbios[@]:((($SLURM_ARRAY_TASK_ID-1) / 40)):1}
idrep=$((($SLURM_ARRAY_TASK_ID-1) % 40 + 1))

julia -p 1 main2.jl 1 $nbio 20 5 1 500 $idrep
