# Notes from getting this thing to build

- Apt-get update wouldn't run but it turned out to be an old version of R which was causing other modules to struggle to update. Fixed this by using a newer version of R.  

- Steps to Convert the docker image to a singularity .sif
1. Push the image to docker hub through VS Code
2. On HPC, use these commands to convert to a .sif
module load containers/apptainer
apptainer build rayshade.sif docker://kmcquil/rayshade:latest

- I ran a test R script that reads data and write out data to make sure everything was working alright. 
- There was one weird error "FATAL:   container creation failed: failed to resolve session directory /localscratch/apptainer/mnt/session: lstat /localscratch/apptainer: no such file or directory". I solved this by just making these directories. Super weird idk

- Here is the example command to run a .sif file 
apptainer exec \
    --pwd /projects/swot/kmcquil/Global_Hillshade \
    --bind /projects/swot/kmcquil/Global_Hillshade \
    --cleanenv \
    /projects/swot/kmcquil/Global_Hillshade/docker/rayshade.sif Rscript src/docker_test.R

--pwd sets the working directory if it is different from the directory where you launch the script 
--bind binds the outside directory with the inside directory. The format is [/source:/dest] but to bind to to the root then you can leave it blank 
--cleanenv is just good practice bc there can be so many wacky variables that we don't want to import into our container
