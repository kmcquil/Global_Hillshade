#!/bin/bash
###########################################################################
## environment & variable setup
####### job customization
#SBATCH --job-name dwnld
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30GB
#SBATCH --time 3:00:00
#SBATCH -p normal_q
#SBATCH -A swot
####### end of job customization
# end of environment & variable setup

module load apptainer/1.4.0
apptainer exec \
    --pwd /projects/swot/kmcquil/Global_Hillshade \
    --bind /projects/swot/kmcquil/Global_Hillshade \
    --cleanenv \
    /projects/swot/kmcquil/Global_Hillshade/docker/rayshade_km.sif python3 src/retile_dem.py