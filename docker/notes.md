# Notes from getting this thing to build

- The rayshade_r_py.sif was built from a Dockerfile using the original Rayshade package (devtools::install_github("tylermorganwall/rayshader")).
- The rayshade_km.sif was built from the Dockerfile in this folder and uses my updated version of the Rayshade package.

### Create a dockerfile 
- An Ubuntu R image is used as the base for the image. 
- Necessary R libraries are installed from requirements.R.

### Steps to create the docker image
- Create the Dockerfile using an Ubuntu R image and install the necessary R libraries from requirements.R
- Build the docker image using the command
```
    docker build -t rayshade .
```
- Push the image to the docker hub. I used VS code. The address is docker://kmcquil/rayshade:latest

Example command to run an rscript from a docker image
```
docker run --mount type=bind,source=C:/Users/kmcquil/Documents/Global_Hillshade,target=/app rayshade Rscript src/docker_test.R
```
Example command to run a python script from a docker image
```
docker run --mount type=bind,source=C:/Users/kmcquil/Documents/Global_Hillshade,target=/app rayshade python3 src/docker_test.py
```

### Steps to create the singularity container (".sif")
- Login to VT cluster and cd to folder where you want to build the .sif file
- Load the necessary modules and build from the docker file. Use the following commands
```
module load containers/apptainer
apptainer build rayshade.sif docker://kmcquil/rayshade:latest
```
- I ran into a weird error ("FATAL:container creation failed: failed to resolve session directory /localscratch/apptainer/mnt/session: lstat /localscratch/apptainer: no such file or directory".). The way to solve this is just by making the directory. Obvious in hindsight. 

Example command to run a .sif file 
```
apptainer exec \
    --pwd /projects/swot/kmcquil/Global_Hillshade \
    --bind /projects/swot/kmcquil/Global_Hillshade \
    --cleanenv \"FATAL:   container creation failed: failed to resolve session directory /localscratch/apptainer/mnt/session: lstat /localscratch/apptainer: no such file or directory".
    /projects/swot/kmcquil/Global_Hillshade/docker/rayshade.sif Rscript src/docker_test.R
```
- pwd sets the working directory if it is different from the directory where you launch the script 
- bind binds the outside directory with the inside directory. The format is [/source:/dest] but to bind to to the root then you can leave it blank 
- cleanenv is just good practice bc there can be so many wacky variables that we don't want to import into our container