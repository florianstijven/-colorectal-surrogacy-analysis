#!/bin/bash
#
#SBATCH --job-name=colorectal
#SBATCH --ntasks=1 --cpus-per-task=36 --mem-per-cpu=5000M
#SBATCH --time=00:45:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A intro_vsc35197

module unuse /apps/leuven/skylake/2018a/modules/all
module use /apps/leuven/skylake/2021a/modules/all

module purge
module load GSL
module load CMake
module load R

Rscript colorectal-model-fit-main.R 36 $SLURM_SUBMIT_DIR



