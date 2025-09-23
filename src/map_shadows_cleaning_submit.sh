#!/bin/bash
###########################################################################
## environment & variable setup
####### job customization
#SBATCH --job-name shadows
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=17GB
#SBATCH --time 100:00:00
#SBATCH -p normal_q
#SBATCH -A swot
#SBATCH --array=1-30
####### end of job customization
# end of environment & variable setup

##module load containers/apptainer
module load apptainer/1.4.0

apptainer exec \
    --pwd /projects/swot/kmcquil/Global_Hillshade \
    --bind /projects/swot/kmcquil/Global_Hillshade \
    --cleanenv \
    /projects/swot/kmcquil/Global_Hillshade/docker/rayshade_km.sif Rscript src/map_shadows.R "data/dem_tiles_cleanup/tiles${SLURM_ARRAY_TASK_ID}.csv" 6 "data/outputs/shadows"