#!/bin/bash
#SBATCH --account=def-ksmccann
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-100
#SBATCH --mail-user=kcazelle@uoguelph.ca
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

Rscript ./launch.R $SLURM_ARRAY_TASK_ID