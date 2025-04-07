# Start from an ubuntu R image version 4.4.1
FROM rocker/r-ver:4.4.1

# Update/install all necessary ubuntu stuff
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libpng-dev libjpeg-dev libfreetype6-dev libglu1-mesa-dev libgl1-mesa-dev pandoc zlib1g-dev libicu-dev libgdal-dev gdal-bin libgeos-dev libproj-dev libmysqlclient-dev libudunits2-dev

# Copy our R script with requirements to the container
COPY requirements.R .

# Install R libraries 
RUN Rscript requirements.R

# Not setting up directories bc they don't seem to be maintained once converted to .sif

# Command to build the image
# docker build -t rayshade .

# Command to run a container from the image, give an rscript as a command, and bind the folders 
# docker run --mount type=bind,source=C:/Users/kmcquil/Documents/Global_Hillshade,target=/ rayshade Rscript src/docker_test.R