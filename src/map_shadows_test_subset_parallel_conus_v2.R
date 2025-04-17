##############################################################################
##############################################################################
# Calculate shadows for a subset of tiles in NA to test the full process
##############################################################################
##############################################################################

# Load packages
library(terra)
library(suntools)
library(data.table)
library(rayshader)
library(sf)
library(foreach)
library(doParallel)

##############################################################################
# Set arguments 
# args[1] is a filepath for a csv with list of tilenames 
##############################################################################
args = commandArgs(trailingOnly=TRUE)
files <- fread(args[1])
n_cores <- as.numeric(args[2]) - 1

##############################################################################
# Functions
##############################################################################

convert_wgs_to_utm <- function(lon, lat){
    # Based on the lat and long, return the best utm epsg-code
    utm_band <- as.character(floor((lon + 180) / 6) + 1)
    utm_band <- sprintf("%02s", utm_band)
    if(lat>0){
        epsg_code <- paste0('326', utm_band)
        return(epsg_code)
    }else{
        epsg_code <- paste0('327', utm_band)
        return(epsg_code)
    }
}


calculate_shadows <- function(
                            dem,
                            lon, 
                            lat, 
                            out_filename,
                            month, 
                            hour, 
                            size, 
                            res_meters){

    # Convert dem, lon, and lat to matrices in the ray shader format
    # Because rayshader reverses the rows of the dem before calculating the shadows
    # I also need to reverse the rows of these so that the solar position and elevation
    # are accurate for the chunk that is being processed
    fix <- function(r){
        r_mat <- terra::as.matrix(r, wide=TRUE)
        r_mat <- r_mat[,ncol(r_mat):1]
        r <- rast(vals=r_mat, crs=crs(dem), extent=ext(dem), resolution=res(dem))
        r_mat <- raster_to_matrix(r)
    }
    
    lon_mat <- fix(lon)
    lat_mat <- fix(lat)
    dem_mat_f <- fix(dem)
    dem_mat <- raster_to_matrix(dem)
    # Create a matrix to store the shadow calcs 
    cache_mask <- matrix(NA, nrow=nrow(dem_mat), ncol=ncol(dem_mat))
   
    # Create the subsets 
    i_seq <- seq(1, nrow(dem_mat), size)
    j_seq <- seq(1, ncol(dem_mat), size)

    # Loop through chunks of cells, calculate the solar position, and then ray trace the shadows
    for(i in i_seq){
        for(j in j_seq){
            # Set up chunks
            i_end = i+size-1
            if(i_end>nrow(dem_mat)){i_end=nrow(dem_mat)}
            j_end = j+size-1
            if(j_end>ncol(dem_mat)){j_end=ncol(dem_mat)}
            
            # Calculate mean elevation in the subset for solar position calculation
            elevation <- mean(dem_mat_f[i:i_end,j:j_end], na.rm=TRUE)

            # Create shadow mask. Only raytrace shadows in cells = 1, skip cells=0
            shadow_mask <- matrix(0, nrow=nrow(dem_mat),ncol=ncol(dem_mat))
            shadow_mask[i:i_end,j:j_end] <- 1
            
            # Calculate solar position at that location
            lon_p <- mean(lon_mat[i:i_end,j:j_end], na.rm=TRUE)
            lat_p <- mean(lat_mat[i:i_end,j:j_end], na.rm=TRUE)
            solar_pos <- solarpos(matrix(c(lon_p, lat_p), nrow = 1),
                    as.POSIXct(paste0("2023-", month, "-01 ", hour, ":00:00"), tz = "EST"),
                    elev = elevation)
            solar_azimuth = as.numeric(solar_pos[1,1])
            solar_elevation = as.numeric(solar_pos[1,2])

            # Perform ray tracing at the specified location
            shadows = ray_shade(
                dem_mat,
                sunaltitude = solar_elevation,
                sunangle = solar_azimuth,
                maxsearch = NULL,
                zscale = res_meters,
                multicore = FALSE,
                cache_mask = NULL,
                shadow_cache = shadow_mask,
                progbar = FALSE,
                anglebreaks = NULL
            )
            
            cache_mask[i:i_end,j:j_end] <- shadows[i:i_end,j:j_end]    
        }
    }

    # Put in the right format
    shadows = t(apply(cache_mask, 2, rev))
    shadow_rast <- rast(vals=shadows, crs=crs(dem), extent=ext(dem), resolution=res(dem))
    writeRaster(shadow_rast, out_filename, overwrite=TRUE)

    return(shadow_rast)
}

apply_calculate_shadows <- function(tif_filename, month, hour, size, out){
    # Load the raster
    dem <- rast(tif_filename)

    # Create matching rasters of lattidue and longitude
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')

    # Project the raster to utm to get the resolution in meters
    utm_crs <- paste0(
        "EPSG:", 
        convert_wgs_to_utm(mean(values(lon)), mean(values(lat)))
        )
    dem_projected <- project(dem, utm_crs) 
    res_meters <- res(dem_projected)[1] 

    # Set the name to write the tif
    out_filename <- paste0(
        out,
        as.character(size), 
        "_", 
        month, 
        "_",
        hour,
        "_",
        basename(tif_filename)
        )

    # Calcualte the shadows
    calculate_shadows(
        dem,
        lon, 
        lat, 
        out_filename,
        month, 
        hour, 
        size, 
        res_meters)

    return()

    }


##############################################################################
# Map shadows
##############################################################################

##############################################################################
# Use four tiles with 2000x2000 chunks and time how long each takes 
# to calculate all of the time steps
##############################################################################
# Set up parallel
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Set directory to store the results
dir.create("data/outputs/test_subset_2000_parallel_v2")
# Create a dt of month/hour combos
months <- c(7)
hours <-c(10, 15)
dts <- data.table(expand.grid(months, hours))
colnames(dts) <- c("month", "hour")
# Set the size of chunks to use for the solar position
size <- 2000
# Set the path to write out shadow maps
out <- "data/outputs/test_subset_2000_parallel_v2/shadow_"

# Calculate shadows and time it
start_time <- Sys.time()

foreach(i=1:nrow(files)) %:%
    foreach(j=1:nrow(dts)) %dopar% {
    #foreach(j=1:5) %dopar% {
        
        # Load packages
        library(terra)
        library(suntools)
        library(data.table)
        library(rayshader)
        library(sf)

        file <- files$file[i]
        print(file)
        apply_calculate_shadows(file, as.character(dts$month[j]), as.character(dts$hour[j]), size, out)

    }

end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

stopCluster(cl)