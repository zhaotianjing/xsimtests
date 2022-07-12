#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500
#SBATCH --time=90:00:00
#SBATCH --partition=high
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu



module load julia

julia save_var.jl