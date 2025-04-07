#!/bin/bash
###########################################################################
## environment & variable setup
####### job customization
#SBATCH --job-name dwnld
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --time 3:00:00
#SBATCH -p normal_q
#SBATCH -A swot
####### end of job customization
# end of environment & variable setup

module load site/tinkercliffs-rome/easybuild/setup
module load site/tinkercliffs/easybuild/setup
module load Miniconda3/23.10.0-1
source activate /home/kmcquil/env/swot_normalq

# Run my script
python download_data.py