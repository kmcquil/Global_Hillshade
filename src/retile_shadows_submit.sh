#!/bin/bash
###########################################################################
## environment & variable setup
####### job customization
#SBATCH --job-name retile
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --time 5:00:00
#SBATCH -p normal_q
#SBATCH -A swot
#SBATCH --array=1-30
####### end of job customization
# end of environment & variable setup

module load apptainer/1.4.0
apptainer exec \
    --pwd /projects/swot/kmcquil/Global_Hillshade \
    --bind /projects/swot/kmcquil/Global_Hillshade \
    --cleanenv \
    /projects/swot/kmcquil/Global_Hillshade/docker/rayshade_km.sif python3 src/retile_shadows.py "data/shadow_tiles/tiles${SLURM_ARRAY_TASK_ID}.csv"
