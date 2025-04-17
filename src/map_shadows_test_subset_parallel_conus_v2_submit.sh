#!/bin/bash
###########################################################################
## environment & variable setup
####### job customization
#SBATCH --job-name shadows
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=12GB
#SBATCH --time 10:00:00
#SBATCH -p normal_q
#SBATCH -A swot
#SBATCH --array=1-10
####### end of job customization
# end of environment & variable setup

module load containers/apptainer
apptainer exec \
    --pwd /projects/swot/kmcquil/Global_Hillshade \
    --bind /projects/swot/kmcquil/Global_Hillshade \
    --cleanenv \
    /projects/swot/kmcquil/Global_Hillshade/docker/rayshade_km.sif Rscript src/map_shadows_test_subset_parallel_conus_v2.R "data/tiles_conus/tiles${SLURM_ARRAY_TASK_ID}.csv" 6