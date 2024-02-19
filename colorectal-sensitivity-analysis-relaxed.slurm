#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=36
#SBATCH --time=01:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1
module list

Rscript -e "renv::restore()"
Rscript colorectal-sensitivity-analysis-relaxed.R 36



