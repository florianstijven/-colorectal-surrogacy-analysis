#!/bin/bash
#
#SBATCH --job-name=colorectal
#SBATCH --ntasks=1 --cpus-per-task=36 --mem-per-cpu=5000M
#SBATCH --time=02:30:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A intro_vsc35197

module unuse /apps/leuven/skylake/2018a/modules/all
module use /apps/leuven/skylake/2021a/modules/all

module purge
module load GSL
module load CMake
module load R

Rscript -e 'renv::run(script = "colorectal-sensitivity-analysis-relaxed.R", project = "/data/leuven/351/vsc35197/Surrogacy/-colorectal-surrogacy-analysis")' 




