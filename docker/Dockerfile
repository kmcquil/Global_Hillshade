# Start from an ubuntu R image version 4.4.1
FROM rocker/r-ver:4.4.1

# Set the working directory
WORKDIR /app

# Update/install all necessary ubuntu stuff
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libglu1-mesa-dev \
    libgl1-mesa-dev \
    pandoc zlib1g-dev \
    libicu-dev \
    libgdal-dev \
    gdal-bin \
    libgeos-dev \
    libproj-dev \
    libmysqlclient-dev \
    libudunits2-dev

# Copy our R script with requirements to the container
COPY requirements.R .

# Install R libraries 
RUN Rscript requirements.R

# Install Python too 
RUN apt-get install -y --no-install-recommends \
    python3 \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip
COPY requirements.txt .
RUN pip install -r requirements.txt
