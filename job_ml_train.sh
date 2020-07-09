#!/bin/bash
#SBATCH --account=def-ksmccann
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-100
#SBATCH --mail-user=kcazelle@uoguelph.ca
#SBATCH --mail-type=FAIL

mkdir res_ml
julia -p 8 src/main_train.jl