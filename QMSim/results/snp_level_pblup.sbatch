#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500
#SBATCH --time=9-99:00:00
#SBATCH --partition=high
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu



module load julia

srun julia snp_level_pblup.jl $1