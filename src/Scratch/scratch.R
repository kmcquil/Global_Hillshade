############################################################################
# I need to figure out what's happening with some of this weird looking 
# mosaicing. 
# flordia looks weird
############################################################################

# Florida weirdness
# n25w085
# n30w085

# Weirdness on east coast a little north of florida
# n30w080
# n30w085

# Try calculating the shadows without an elevation 

# Load packages
library(terra)
library(suntools)
library(data.table)
library(rayshader)
library(sf)
library(foreach)
library(doParallel)

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
                    as.POSIXct(paste0("2023-", month, "-01 ", hour, ":00:00"), tz = "EST"))
                    #elev = elevation)
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


# Florida weirdness
# n25w085
# n30w085

# Weirdness on east coast a little north of florida
# n30w080
# n30w085

setwd("C:/Users/kmcquil/Documents/Global_Hillshade")

files <- data.table("file"=c("data/raw/merit_retile/n25w085_elv.tif", "data/raw/merit_retile/n30w085_elv.tif", "data/raw/merit_retile/n30w080_elv.tif"))

months <- c(7)
hours <-c(15)
dts <- data.table(expand.grid(months, hours))
colnames(dts) <- c("month", "hour")
# Set the size of chunks to use for the solar position
size <- 2000
# Set the path to write out shadow maps
out <- "data/outputs/scratch/shadow_"

cl <- makeCluster(3)
registerDoParallel(cl)

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
stopCluster(cl)



##########################################################################
# Test the solar position 


bad = rast("data/raw/merit_retile/n30w080_elv.tif")
plot(bad)

good = rast("data/raw/merit_retile/n30w085_elv.tif")
plot(good)


##########################################################################
# Lets take my weirdness out of it. If i give them both the 
# same solar position for the full tile, will it match up
# if it does, then im doing some jenky shit with my code with subset


reg_shadows <- function(file, outname){
    dem <- rast(file)
    dem_mat <- raster_to_matrix(dem)
    shadows <- ray_shade(
                dem_mat,
                sunaltitude = 54,
                sunangle = 262,
                maxsearch = NULL,
                zscale = 90,
                multicore = FALSE,
                cache_mask = NULL,
                shadow_cache = NULL,
                progbar = FALSE,
                anglebreaks = NULL
            )
    shadows_rast = t(apply(shadows, 2, rev))
    shadows_rast <- rast(vals=shadows_rast, crs=crs(dem), extent=ext(dem), resolution=res(dem))
    writeRaster(shadows_rast, paste0("data/outputs/scratch/", outname), overwrite=TRUE)

    return(shadows)
}


bad_shadows <- reg_shadows("data/raw/merit_retile/n30w080_elv.tif", "test_bad.tif")
good_shadows <- reg_shadows("data/raw/merit_retile/n30w085_elv.tif", "test_good.tif")

solar_pos <- solarpos(matrix(c(-80, 34), nrow = 1),
                    as.POSIXct(paste0("2023-", "07", "-01 ", "15", ":00:00"), tz = "EST"))
                    #elev = elevation)
solar_azimuth = as.numeric(solar_pos[1,1])
solar_elevation = as.numeric(solar_pos[1,2])




##########################################################################
# Now try two small subsets that don't have any nas in the right corner and see 
library(sf)
reg_shadows <- function(file, outname, poly){
    dem <- rast(file)
    poly <- st_read(poly)
    dem <- crop(dem, poly)
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')
    lon_mat <- raster_to_matrix(lon) 
    lat_mat <- raster_to_matrix(lat)
    lon_p <- mean(lon_mat, na.rm=TRUE)
    lat_p <- mean(lat_mat, na.rm=TRUE)
    solar_pos <- solarpos(matrix(c(lon_p, lat_p), nrow = 1),
                    as.POSIXct(paste0("2023-", "07", "-01 ", "15", ":00:00"), tz = "EST"))
                    #elev = elevation)
    solar_azimuth = as.numeric(solar_pos[1,1])
    solar_elevation = as.numeric(solar_pos[1,2])
    print(outname)
    print(paste0("solar_azimuth=", solar_azimuth))
    print(paste0("solar_elevation=", solar_elevation))
    dem_mat <- raster_to_matrix(dem)
    shadows <- ray_shade(
                dem_mat,
                sunaltitude = solar_elevation,
                sunangle = solar_azimuth,
                maxsearch = NULL,
                zscale = 90,
                multicore = FALSE,
                cache_mask = NULL,
                shadow_cache = NULL,
                progbar = FALSE,
                anglebreaks = NULL
            )
    shadows_rast = t(apply(shadows, 2, rev))
    shadows_rast <- rast(vals=shadows_rast, crs=crs(dem), extent=ext(dem), resolution=res(dem))
    writeRaster(shadows_rast, paste0("data/outputs/scratch/", outname), overwrite=TRUE)

    return(shadows)
}

bad_shadows <- reg_shadows("data/raw/merit_retile/n30w080_elv.tif", 
                            "test_bad.tif", 
                            "data/outputs/scratch/poly_test_bad.shp")
good_shadows <- reg_shadows("data/raw/merit_retile/n30w085_elv.tif", 
                            "test_good.tif", 
                            "data/outputs/scratch/poly_test_good.shp")









##########################################################################
# Now try different angle breaks with a smaller total angle range to 
# see if that helps



library(sf)
reg_shadows <- function(file, outname, poly){
    dem <- rast(file)
    poly <- st_read(poly)
    dem <- crop(dem, poly)
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')
    lon_mat <- raster_to_matrix(lon) 
    lat_mat <- raster_to_matrix(lat)
    lon_p <- mean(lon_mat, na.rm=TRUE)
    lat_p <- mean(lat_mat, na.rm=TRUE)
    solar_pos <- solarpos(matrix(c(lon_p, lat_p), nrow = 1),
                    as.POSIXct(paste0("2023-", "07", "-01 ", "15", ":00:00"), tz = "EST"))
                    #elev = elevation)
    solar_azimuth = as.numeric(solar_pos[1,1])
    solar_elevation = as.numeric(solar_pos[1,2])
    angle_of_sun = 0.01 # it should actually realistically be 0.533
    anglebreaks = seq(max(0,solar_elevation-angle_of_sun/2), min(90,solar_elevation+angle_of_sun/2), length.out = 1)
    print(outname)
    print(paste0("solar_azimuth=", solar_azimuth))
    print(paste0("solar_elevation=", solar_elevation))
    dem_mat <- raster_to_matrix(dem)
    shadows <- ray_shade(
                dem_mat,
                #sunaltitude = solar_elevation,
                sunangle = solar_azimuth,
                maxsearch = NULL,
                zscale = 90,
                multicore = FALSE,
                cache_mask = NULL,
                shadow_cache = NULL,
                progbar = FALSE,
                anglebreaks = anglebreaks
            )
    shadows_rast = t(apply(shadows, 2, rev))
    shadows_rast <- rast(vals=shadows_rast, crs=crs(dem), extent=ext(dem), resolution=res(dem))
    writeRaster(shadows_rast, paste0("data/outputs/scratch/", outname), overwrite=TRUE)

    return(shadows)
}

bad_shadows <- reg_shadows("data/raw/merit_retile/n30w080_elv.tif", 
                            "test_bad.tif", 
                            "data/outputs/scratch/poly_test_bad.shp")
good_shadows <- reg_shadows("data/raw/merit_retile/n30w085_elv.tif", 
                            "test_good.tif", 
                            "data/outputs/scratch/poly_test_good.shp")




##########################################################################
# Now I want to calculate the shadows with them combined
# Merge those elevation tiles together and then calculate the shadow

r1 <- rast("data/raw/merit_retile/n30w080_elv.tif")
r2 <- rast("data/raw/merit_retile/n30w085_elv.tif")

s <- sprc(r1, r2)
m <- merge(s)
writeRaster(m, "data/outputs/scratch/combo_rast.tif")

reg_shadows <- function(file, outname, poly){
    dem <- rast(file)
    poly <- st_read(poly)
    dem <- crop(dem, poly)
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')
    lon_mat <- raster_to_matrix(lon) 
    lat_mat <- raster_to_matrix(lat)
    lon_p <- mean(lon_mat, na.rm=TRUE)
    lat_p <- mean(lat_mat, na.rm=TRUE)
    solar_pos <- solarpos(matrix(c(lon_p, lat_p), nrow = 1),
                    as.POSIXct(paste0("2023-", "07", "-01 ", "15", ":00:00"), tz = "EST"))
                    #elev = elevation)
    solar_azimuth = as.numeric(solar_pos[1,1])
    solar_elevation = as.numeric(solar_pos[1,2])
    print(outname)
    print(paste0("solar_azimuth=", solar_azimuth))
    print(paste0("solar_elevation=", solar_elevation))
    dem_mat <- raster_to_matrix(dem)
    shadows <- ray_shade(
                dem_mat,
                sunaltitude = solar_elevation,
                sunangle = solar_azimuth,
                maxsearch = NULL,
                zscale = 90,
                multicore = FALSE,
                cache_mask = NULL,
                shadow_cache = NULL,
                progbar = FALSE,
                anglebreaks = NULL
            )
    shadows_rast = t(apply(shadows, 2, rev))
    shadows_rast <- rast(vals=shadows_rast, crs=crs(dem), extent=ext(dem), resolution=res(dem))
    writeRaster(shadows_rast, paste0("data/outputs/scratch/", outname), overwrite=TRUE)

    return()
}

reg_shadows("data/outputs/scratch/combo_rast.tif", 
                            "test_combo.tif", 
                            "data/outputs/scratch/poly_test_combo.shp")




##########################################################################
# Now I want to try setting lambert to FALSE since I think that 
# is attempting to account for the general light intensity 

reg_shadows <- function(file, outname, poly){
    dem <- rast(file)
    poly <- st_read(poly)
    dem <- crop(dem, poly)
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')
    lon_mat <- raster_to_matrix(lon) 
    lat_mat <- raster_to_matrix(lat)
    lon_p <- mean(lon_mat, na.rm=TRUE)
    lat_p <- mean(lat_mat, na.rm=TRUE)
    solar_pos <- solarpos(matrix(c(lon_p, lat_p), nrow = 1),
                    as.POSIXct(paste0("2023-", "07", "-01 ", "15", ":00:00"), tz = "EST"))
                    #elev = elevation)
    solar_azimuth = as.numeric(solar_pos[1,1])
    solar_elevation = as.numeric(solar_pos[1,2])
    print(outname)
    print(paste0("solar_azimuth=", solar_azimuth))
    print(paste0("solar_elevation=", solar_elevation))
    dem_mat <- raster_to_matrix(dem)
    shadows <- ray_shade(
                dem_mat,
                sunaltitude = solar_elevation,
                sunangle = solar_azimuth,
                maxsearch = NULL,
                zscale = 90,
                multicore = FALSE,
                cache_mask = NULL,
                shadow_cache = NULL,
                progbar = FALSE,
                anglebreaks = NULL,
                lambert=FALSE
            )
    shadows_rast = t(apply(shadows, 2, rev))
    shadows_rast <- rast(vals=shadows_rast, crs=crs(dem), extent=ext(dem), resolution=res(dem))
    writeRaster(shadows_rast, paste0("data/outputs/scratch/", outname), overwrite=TRUE)

    return()
}

reg_shadows("data/outputs/scratch/combo_rast.tif", 
                            "test_combo_lambert_false.tif", 
                            "data/outputs/scratch/poly_test_combo.shp")




##########################################################################
# Last thing is to calculate shadows after combing both tiles 
# without cropping to show how it should look


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
                    as.POSIXct(paste0("2023-", month, "-01 ", hour, ":00:00"), tz = "EST"))
                    #elev = elevation)
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

setwd("C:/Users/kmcquil/Documents/Global_Hillshade")

# Set the path to write out shadow maps
out <- "data/outputs/scratch/shadow_"

file <- "data/outputs/scratch/combo_rast.tif"
apply_calculate_shadows(file, "07", "15", 2000, out)






#################################################################################
#################################################################################
# The lambda rescaling is the issue
# Figure out how to scale in a way that is smarter 
#################################################################################
#################################################################################

library(devtools)
devtools::install_github("kmcquil/rayshader_km")

library(terra)
library(suntools)
library(data.table)
library(rayshader)
library(sf)

file1 <- "data/raw/merit_retile/n30w085_elv.tif" # this is the full one 
dem <- rast(file1)
lon <- init(dem, 'x')
lat <- init(dem, 'y')
lon_mat <- raster_to_matrix(lon) 
lat_mat <- raster_to_matrix(lat)
lon_p <- mean(lon_mat, na.rm=TRUE)
lat_p <- mean(lat_mat, na.rm=TRUE)
solar_pos <- solarpos(matrix(c(lon_p, lat_p), nrow = 1),
                    as.POSIXct(paste0("2023-", "07", "-01 ", "15", ":00:00"), tz = "EST"))
                    #elev = elevation)
solar_azimuth = as.numeric(solar_pos[1,1])
solar_elevation = as.numeric(solar_pos[1,2])
print(paste0("solar_azimuth=", solar_azimuth))
print(paste0("solar_elevation=", solar_elevation))
dem_mat <- raster_to_matrix(dem)

# Ray shade code to get ready for lamb shade
heightmap <- dem_mat
sunaltitude=solar_elevation
sunangle=solar_azimuth
maxsearch=NULL
lambert=TRUE
zscale=90
multicore = FALSE
cache_mask = NULL
shadow_cache=NULL
progbar=interactive() 
anglebreaks = NULL

if(is.null(anglebreaks)) {
    anglebreaks = seq(max(0,sunaltitude-0.533/2), min(90,sunaltitude+0.533/2), length.out = 10)
}
if(all(anglebreaks <= 0)) {
    return(matrix(0,nrow=nrow(heightmap),ncol=ncol(heightmap)))
}
if(is.null(maxsearch)) {
    maxsearch = (max(heightmap,na.rm=TRUE) - min(heightmap,na.rm=TRUE))/(zscale*sinpi(min(anglebreaks[anglebreaks > 0])/180))
}
anglebreaks = anglebreaks[order(anglebreaks)]
anglebreaks_rad = anglebreaks*pi/180
sunangle_rad = pi+sunangle*pi/180
originalheightmap = heightmap

# Start of lambshade
test1 = lamb_shade(originalheightmap, sunaltitude = mean(anglebreaks), sunangle = sunangle, zscale = 90)
min(test1, na.rm=TRUE)
max(test1, na.rm=TRUE)


file2 <- "data/raw/merit_retile/n30w080_elv.tif" # this is the small one 
dem <- rast(file2)
lon <- init(dem, 'x')
lat <- init(dem, 'y')
lon_mat <- raster_to_matrix(lon) 
lat_mat <- raster_to_matrix(lat)
lon_p <- mean(lon_mat, na.rm=TRUE)
lat_p <- mean(lat_mat, na.rm=TRUE)
solar_pos <- solarpos(matrix(c(lon_p, lat_p), nrow = 1),
                    as.POSIXct(paste0("2023-", "07", "-01 ", "15", ":00:00"), tz = "EST"))
                    #elev = elevation)
solar_azimuth = as.numeric(solar_pos[1,1])
solar_elevation = as.numeric(solar_pos[1,2])
print(paste0("solar_azimuth=", solar_azimuth))
print(paste0("solar_elevation=", solar_elevation))
dem_mat <- raster_to_matrix(dem)

# Ray shade code to get ready for lamb shade
heightmap <- dem_mat
sunaltitude=solar_elevation
sunangle=solar_azimuth
maxsearch=NULL
lambert=TRUE
zscale=90
multicore = FALSE
cache_mask = NULL
shadow_cache=NULL
progbar=interactive() 
anglebreaks = NULL

if(is.null(anglebreaks)) {
    anglebreaks = seq(max(0,sunaltitude-0.533/2), min(90,sunaltitude+0.533/2), length.out = 10)
}
if(all(anglebreaks <= 0)) {
    return(matrix(0,nrow=nrow(heightmap),ncol=ncol(heightmap)))
}
if(is.null(maxsearch)) {
    maxsearch = (max(heightmap,na.rm=TRUE) - min(heightmap,na.rm=TRUE))/(zscale*sinpi(min(anglebreaks[anglebreaks > 0])/180))
}
anglebreaks = anglebreaks[order(anglebreaks)]
anglebreaks_rad = anglebreaks*pi/180
sunangle_rad = pi+sunangle*pi/180
originalheightmap = heightmap

# Start of lambshade
test2 = lamb_shade(originalheightmap, sunaltitude = mean(anglebreaks), sunangle = sunangle, zscale = 90)
min(test2, na.rm=TRUE)
max(test2, na.rm=TRUE)










































library(devtools)
devtools::install_github("kmcquil/rayshader_km")

library(terra)
library(suntools)
library(data.table)
library(rayshader)
library(sf)


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
                            outdir,
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
    cache_mask_rt <- matrix(NA, nrow=nrow(dem_mat), ncol=ncol(dem_mat))
    cache_mask_ls <- matrix(NA, nrow=nrow(dem_mat), ncol=ncol(dem_mat))
    cache_mask_sh <- matrix(NA, nrow=nrow(dem_mat), ncol=ncol(dem_mat))

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
            
            sh = shadows[[1]]
            rt = shadows[[2]]
            ls = shadows[[3]]
            cache_mask_rt[i:i_end,j:j_end] <- rt[i:i_end,j:j_end]
            cache_mask_ls[i:i_end,j:j_end] <- ls[i:i_end,j:j_end]
            cache_mask_sh[i:i_end,j:j_end] <- sh[i:i_end,j:j_end]
              
        }
    }
    # Function to write out
    transform_and_write <- function(cache_mask, suf){
        shadows = t(apply(cache_mask, 2, rev))
        shadow_rast <- rast(vals=shadows, crs=crs(dem), extent=ext(dem), resolution=res(dem))
        shadow_rast <- mask(shadow_rast, dem)
        # Replace NA values with 2.55 This will be converted to 255 when we multiply by 100 
        shadow_rast <- subst(shadow_rast, NA, 2.55)
        shadow_rast <- round(shadow_rast*100, 0)
        out_fp <- paste0(outdir, "/", suf, "/", out_filename)
        writeRaster(shadow_rast, out_fp, overwrite=TRUE, datatype="INT8U", gdal=c("COMPRESS=DEFLATE"), NAflag=255)
    }
    transform_and_write(cache_mask_rt, "raytrace")
    transform_and_write(cache_mask_ls, "lambshade")
    transform_and_write(cache_mask_sh, "combo")
}

apply_calculate_shadows <- function(tif_filename, month, hour, size, outdir){
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
        "shadow_",
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
        outdir,
        out_filename,
        month, 
        hour, 
        size, 
        res_meters)

    return()

    }


setwd("C:/Users/kmcquil/Documents/Global_Hillshade")

outdir <- "data/outputs/scratch"
# Make sure directory to store results exists, and make it if it does not
dir.create(outdir)
# Make sure the subdirectories exist 
dir.create(paste0(outdir, "/raytrace"))
dir.create(paste0(outdir, "/lambshade"))
dir.create(paste0(outdir, "/combo"))

# Should add code to make sure the sub folders exist for the out
# Set the path to write out shadow maps
file1 <- "data/raw/merit_retile/n30w085_elv.tif" # this is the full one 

apply_calculate_shadows(file1, "07", "8", 3000, outdir)
