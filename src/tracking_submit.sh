#!/bin/bash
###########################################################################
## environment & variable setup
####### job customization
#SBATCH --job-name cleaning
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --time 00:10:00
#SBATCH -p normal_q
#SBATCH -A swot
####### end of job customization
# end of environment & variable setup

module load R
Rscript tracking.R
