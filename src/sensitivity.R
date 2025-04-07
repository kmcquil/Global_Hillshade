##############################################################################
##############################################################################
# Calculate shadows for a subset of tiles in NA to test stuff out
##############################################################################
##############################################################################

# Load packages
library(terra)
library(suntools)
library(data.table)
library(terra)
library(rayshader)
library(sf)
library(ggplot2)

home <- "C:/Users/kmcquil/Documents/Global_Hillshade"

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
                            dem_mat,
                            lon_mat, 
                            lat_mat, 
                            out_filename,
                            month, 
                            hour, 
                            size, 
                            res_meters,
                            lambert){

    # Create a matrix to store the shadow calcs 
    cache_mask <- matrix(NA, nrow(dem_mat), ncol(dem_mat))
   
    # Create the subsets 
    i_seq <- seq(1, nrow(dem_mat), size)
    j_seq <- seq(1, ncol(dem_mat), size)
    print(i_seq)
    print(j_seq)

    # Loop through chunks of cells, calculate the solar position, and then ray trace the shadows
    for(i in i_seq){
        for(j in j_seq){
            # Track progress
            print(paste0("i=", i))
            print(paste0("j=", j))

            # Set up chunks
            i_end = i+size-1
            if(i_end>nrow(dem_mat)){i_end=nrow(dem_mat)}
            j_end = j+size-1
            if(j_end>ncol(dem_mat)){j_end=ncol(dem_mat)}
            
            # Check if the subset is empty. If it is, move onto the next chunk
            elevation <- mean(dem_mat[i:i_end,j:j_end], na.rm=TRUE)
            if(is.na(elevation)){next}

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
            print(solar_azimuth)
            print(solar_elevation)
            print(elevation)

            # Perform ray tracing at the specified location
            shadows = ray_shade(
                dem_mat,
                sunaltitude = solar_elevation,
                sunangle = solar_azimuth,
                maxsearch = NULL,
                lambert = lambert,
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

apply_calculate_shadows <- function(tif_filename, month, hour, size, lambert, out="data/outputs/shadow_test_lambert_"){
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

    # Convert dem, lon, and lat to matrices
    dem_mat <- raster_to_matrix(dem)
    lon_mat <- raster_to_matrix(lon)
    lat_mat <- raster_to_matrix(lat)

    # Set the name to write the tif
    if(lambert==TRUE){
        out_filename <- paste0(
            home, 
            "/",
            out,
            as.character(size), 
            "_", 
            month, 
            "_",
            hour,
            "_",
            basename(tif_filename)
            )
    }
    if(lambert==FALSE){
        # Set the name to write out shadow map
        out_filename <- paste0(
            home, 
            "/",
            out, 
            as.character(size), 
            "_", 
            month, 
            "_",
            hour,
            "_",
            basename(tif_filename)
            )
    }

    # Calcualte the shadows
    calculate_shadows(
        dem,
        dem_mat,
        lon_mat, 
        lat_mat, 
        out_filename,
        month, 
        hour, 
        size, 
        res_meters,
        lambert
        )

    return()

    }


# Function to calc solar position 
calculate_sp <- function(r, month, hour){
    # Calculate the solar position based on location and time
    # r: (vector) r[1] = longitude, r[2] = latitude, r[3] = elevation
    # month: (character, len=2)
    # hour: (character, len=2)
    # return vector with solar azimuth, solar elevation
    datetime = as.POSIXct(paste0("2023-", month, "-01 ", hour, ":00:00"), tz = "EST")
    solar_pos <- solarpos(matrix(c(r[1], r[2]), nrow = 1),
        as.POSIXct(paste0("2023-", month, "-01 ", hour, ":00:00"), tz = "EST"), 
        )#elev = r[3])
    solar_azimuth = as.numeric(solar_pos[1,1])
    solar_elevation = as.numeric(solar_pos[1,2])
    return(data.table(datetime = datetime, solar_azimuth=solar_azimuth, solar_elevation=solar_elevation))
}

calculate_sensitivity <- function(tif_filename, date_times, sizes){
    # Load the DEM
    dem <- rast(tif_filename)
    
    # Create lat/long rasters and stack
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')
    s <- c(lon, lat, dem)

    # Calcualte sp in four corners of the subset for each subset size and datetime
    sp_dt = list()
    i = 0
    for(num_pixels in sizes){
        for(d in date_times){
            i = i + 1
            solar_position_ul = calculate_sp(s[[1:3]][1,1], d[1], d[2])
            solar_position_ll = calculate_sp(s[[1:3]][num_pixels,1], d[1], d[2])
            solar_position_ur = calculate_sp(s[[1:3]][1,num_pixels], d[1], d[2])
            solar_position_lr = calculate_sp(s[[1:3]][num_pixels, num_pixels], d[1], d[2])
            dt = rbind(solar_position_ul, solar_position_ll, solar_position_ur, solar_position_lr)
            dt$size = num_pixels
            sp_dt[[i]] = dt
        }
    }
    sp_dt = rbindlist(sp_dt)
    sp_dt$filename <- tif_filename
    return(sp_dt)
}

##############################################################################
# Map shadows
##############################################################################

##############################################################################
# Test using 3000 x 3000 size subsets 
##############################################################################
month <- "07"
hour <- "15"
size <- 3000

files <- list.files(paste0(home, "/data/raw/merit_retile"), pattern=".tif", full.names=TRUE)

start_time <- Sys.time()

for(file in files){
    print(file)
    apply_calculate_shadows(file, month, hour, size)
}

end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

##############################################################################
# Test using 2000 x 2000 size subsets 
##############################################################################

month <- "07"
hour <- "15"
size <- 2000

files <- list.files(paste0(home, "/data/raw/merit_retile"), pattern=".tif", full.names=TRUE)
files_subset <- files[c(c(3,4), c(9,10))]

start_time <- Sys.time()
for(file in files_subset){
    print(file)
    apply_calculate_shadows(file, month, hour, size)
}
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken


##############################################################################
# Test with lambert = FALSE to get binary of shade or no shade
##############################################################################

month <- "07"
hour <- "15"
sizes <-c(2000, 3000)

files <- list.files(paste0(home, "/data/raw/merit_retile"), pattern=".tif", full.names=TRUE)
files_subset <- files[c(c(3,4), c(9,10))]

for(size in sizes){
    for(file in files_subset){
        apply_calculate_shadows(file, month, hour, size, lambert=FALSE)
    }
}


##############################################################################
# Test at different latitudes and different times of day 
##############################################################################

# Two possible size chunks
sizes <-c(2000, 3000)

# three files that are in EST and span a latitude gradient
files_subset <- c(
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n00w075_elv.tif",
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n35w085_elv.tif", 
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n60w075_elv.tif"
    )

for(size in sizes){
    for(file in files_subset){
        apply_calculate_shadows(file, "07", "08", size, lambert=TRUE) # 8 am 
        apply_calculate_shadows(file, "07", "13", size, lambert=TRUE) # 1 pm 
        apply_calculate_shadows(file, "07", "18", size, lambert=TRUE) # 6 pm
    }
} 






##############################################################################
# This is the big sensitivity test that actually matters
# Choose three tiles in the EST time zone at latitudes = 0, 30, 60
# Calculate shadows in those tiles with subset sizes of 1000, 2000, 3000
# at times 9am, 12pm, 3pm, 6pm
# Also, for the subset in the top left corner, calculate the solar position 
# in each corner and save it for that latitude, time, and subset sie
# After we have calcualted all of this, perform a pixel size 
# sensitivity analysis of the shadows and the solar position
# For example, how different is the solar position in degrees based on those factors
# For example, what is the average difference in shadow intensity based on those factors 
##############################################################################

##############################################################################
# Calculate sp and shadows for the given dates, subset sizes, and latitudes
##############################################################################

date_times <- list(c("07", "09"), c("07", "12"), c("07", "15"), c("07", "18"))
sizes <- c(1000, 2000, 3000)
files <- c(
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n00w075_elv.tif",
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n35w085_elv.tif", 
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n60w075_elv.tif"
    )
out <- "data/outputs/sensitivity_test/shadow_"

sensitivity_df <- rbindlist(lapply(files, calculate_sensitivity, date_times=date_times, sizes=sizes))
fwrite(sensitivity_df, paste0(out, "sensitivity.csv"))

for(date_time in date_times){
    for(size in sizes){
        for(file in files){
            apply_calculate_shadows(file, date_time[1], date_time[2], size, lambert=TRUE, out=out) 
        }
    } 
}

##############################################################################
# Assess sensitivity of sp
##############################################################################

# Calculate the difference in solar azimuth and solar elevation by latitude, size, and hour of day
sensitivity_df$hour <- as.numeric(format(as.POSIXct(sensitivity_df$datetime), format = "%H"))
sensitivity_df$Latitude <- as.numeric(substr(basename(sensitivity_df$filename), 2, 3))

calc_diff <- function(column_name){
    diff1 = max(column_name) - min(column_name)
    diff2 = 360-diff1
    return(min(c(diff1, diff2)))   
}
sensitivity <- sensitivity_df[, .(diff_solar_azimuth=calc_diff(solar_azimuth), 
                        diff_solar_elevation=calc_diff(solar_elevation)), 
                        by=.(Latitude, hour, size)]

sensitivity_agg <- sensitivity[,.(min_diff_solar_azimuth=min(diff_solar_azimuth),
                                mean_diff_solar_azimuth=mean(diff_solar_azimuth), 
                                max_diff_solar_azimuth=max(diff_solar_azimuth), 
                                min_diff_solar_elevation=min(diff_solar_elevation),
                                mean_diff_solar_elevation=mean(diff_solar_elevation),
                                max_diff_solar_elevation=max(diff_solar_elevation)), by=size]
sensitivity_agg[, c("size", "min_diff_solar_azimuth", "mean_diff_solar_azimuth", "max_diff_solar_azimuth")]
print(sensitivity_agg[, c("size", "min_diff_solar_elevation", "mean_diff_solar_elevation", "max_diff_solar_elevation")])

# Plot max difference between the four corners according to the size
sensitivity$Latitude <- as.factor(sensitivity$Latitude)
ggplot(sensitivity) + 
    geom_point(aes(x=size, y=diff_solar_azimuth, color=Latitude), size=5) +
    xlab("Subset size (NxN pixel subset)") + 
    ylab("Max difference in solar azimuth between four corners (degrees)") + 
    theme_bw() + 
    facet_wrap(~hour)

ggplot(sensitivity) + 
    geom_point(aes(x=size, y=diff_solar_elevation, color=Latitude), size=5) +
    xlab("Subset size (NxN pixel subset)") + 
    ylab("Max difference in solar elevation between four corners (degrees)") + 
    theme_bw() + 
    facet_wrap(~hour)


##############################################################################
# Assess sensitivity of shadows
##############################################################################

# I have a tif for each latitude, hour, and subset
file_dt <- data.table(file=list.files("data/outputs/sensitivity_test", full.names=TRUE, pattern=".tif"))
file_dt$subset <- substr(basename(file_dt$file), 8, 11)
file_dt$hour <- substr(basename(file_dt$file), 16, 17)
file_dt$latitude <- substr(basename(file_dt$file), 20, 21)

# Calculate the difference in shadow intensity 
# For a given latitude: 
#   For each time, calculate the pixel-wise difference in shadow intensity by subset (1000-2000, 1000-3000, 2000-3000)
#   Make a table with quantiles of the differences for each. This would have 12 rows with columns including 
#   Latitude, Time, Subset difference, 2.5%, 25%, 50%, 75%, 97.5%
calculate_quantiles <- function(r, label, latitude, hour){
    q_dt <- as.data.table(as.data.frame(global(r, quantile, probs=c(0, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 1), na.rm=TRUE)))
    #names(q_dt) <- c("q0", "q01", "q025", "q10", "q25", "q50", "q75", "q90", "q975", "q99", "q1")
    names(q_dt) <- c("0.0", "0.01", "0.025", "0.10", "0.25", "0.50", "0.75", "0.90", "0.975", "0.99", "1.0")
    q_dt$label <- label
    q_dt$latitude <- latitude
    q_dt$hour <- hour
    return(q_dt)
}
latitudes <- c("00", "35", "60")
hours <- c("09", "12", "15", "18")
quantile_dt <- list()
q <- 0
for(lat in latitudes){
    for(h in hours){
        r_1000 <- rast(file_dt[latitude==lat&hour==h&subset==1000,]$file)
        r_2000 <- rast(file_dt[latitude==lat&hour==h&subset==2000,]$file)
        r_3000 <- rast(file_dt[latitude==lat&hour==h&subset==3000,]$file)
        diff_1000_2000 <- abs(r_1000-r_2000)
        diff_1000_3000 <- abs(r_1000-r_3000)
        diff_2000_3000 <- abs(r_2000-r_3000)
        q <- q + 1
        quantile_dt[[q]] <- calculate_quantiles(diff_1000_2000, label="1000-2000", lat, h)
        q <- q + 1
        quantile_dt[[q]] <- calculate_quantiles(diff_1000_3000, label="1000-3000", lat, h)
        q <- q + 1
        quantile_dt[[q]] <- calculate_quantiles(diff_2000_3000, label="2000-3000", lat, h)
    }
}
quantile_dt <- rbindlist(quantile_dt)
quantile_dt_long <- melt(quantile_dt, 
                        id.vars=c("label", "latitude", "hour"), 
                        measure.vars=c("0.0", "0.01", "0.025", "0.10", "0.25", "0.50", "0.75", "0.90", "0.975", "0.99", "1.0"), 
                        variable.name = "quantile", 
                        value.name = "value")
quantile_dt_long$quantile <- as.numeric(as.character(quantile_dt_long$quantile))

summary(quantile_dt)


ggplot(quantile_dt_long[label=="1000-2000"]) + 
    geom_point(aes(x=quantile, y=value, color=hour), size=3) + 
    theme_bw() + 
    xlab("Quantile") + 
    ylab("Difference in shade intensity") + 
    facet_wrap(~latitude) + 
    ggtitle("1000-2000")

ggplot(quantile_dt_long[label=="1000-3000"]) + 
    geom_point(aes(x=quantile, y=value, color=hour), size=3) + 
    theme_bw() + 
    xlab("Quantile") + 
    ylab("Difference in shade intensity") + 
    facet_wrap(~latitude) + 
    ggtitle("1000-3000")

ggplot(quantile_dt_long[label=="2000-3000"]) + 
    geom_point(aes(x=quantile, y=value, color=hour), size=3) + 
    theme_bw() + 
    xlab("Quantile") + 
    ylab("Difference in shade intensity") + 
    facet_wrap(~latitude) + 
    ggtitle("2000-3000")






##############################################################################
# now try two new things 
# 1. Quantify the difference in solar position from the center of a chunk 
# where it is the most accurate with the four corners
# 2. In the bottom left corner, calculate how different that assigned solar position 
# is from the cells 
##############################################################################

calculate_sp_corner_diffs <- function(tif_filename, date_times, sizes){
    # Load the DEM
    dem <- rast(tif_filename)
    
    # Create lat/long rasters and stack
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')
    s <- c(lon, lat, dem)

    # Calcualte sp in four corners of the subset for each subset size and datetime
    sp_dt = list()
    i = 0
    for(num_pixels in sizes){
        for(d in date_times){
            i = i + 1
            c = round(num_pixels/2)
            solar_position_c = calculate_sp(s[[1:3]][c,c], d[1], d[2])
            solar_position_ul = calculate_sp(s[[1:3]][1,1], d[1], d[2])
            solar_position_ll = calculate_sp(s[[1:3]][num_pixels,1], d[1], d[2])
            solar_position_ur = calculate_sp(s[[1:3]][1,num_pixels], d[1], d[2])
            solar_position_lr = calculate_sp(s[[1:3]][num_pixels, num_pixels], d[1], d[2])
            dt = rbind(solar_position_ul, solar_position_ll, solar_position_ur, solar_position_lr)
            dt$az_diff <- NA
            for(w in 1:nrow(dt)){
                d1 <- abs(dt$solar_azimuth[w] - solar_position_c$solar_azimuth[1])
                d2 <- 360-d1
                dt$az_diff[w] <- min(d1, d2)
            }
            solar_position_c$max_solar_azimuth_difference <- max(dt$az_diff)
            solar_position_c$max_solar_elevation_difference <- max(abs(dt$solar_elevation - solar_position_c$solar_elevation))
            solar_position_c$size = num_pixels
            solar_position_c <- solar_position_c[,!c("solar_azimuth", "solar_elevation")]
            sp_dt[[i]] = solar_position_c
        }
    }
    sp_dt = rbindlist(sp_dt)
    sp_dt$filename <- tif_filename
    return(sp_dt)
}

date_times <- list(c("07", "09"), c("07", "12"), c("07", "15"), c("07", "18"))
sizes <- c(1000, 2000, 3000)
files <- c(
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n00w075_elv.tif",
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n35w085_elv.tif", 
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n60w075_elv.tif"
    )
out <- "data/outputs/sensitivity_test"

sp4c_df <- rbindlist(lapply(files, calculate_sp_corner_diffs, date_times=date_times, sizes=sizes))
sp4c_df$hour <- as.numeric(format(as.POSIXct(sp4c_df$datetime), format = "%H"))
sp4c_df$Latitude <- as.numeric(substr(basename(sp4c_df$filename), 2, 3))
fwrite(sp4c_df, paste0(out, "/sp4c_df.csv"))

# Plot max difference between the four corners according to the size 
sp4c_df$Latitude <- as.factor(sp4c_df$Latitude)
ggplot(sp4c_df) + 
    geom_point(aes(x=size, y=max_solar_azimuth_difference, color=Latitude), size=5) +
    xlab("Subset size (NxN pixel subset)") + 
    ylab("Max difference in solar azimuth between center and four corners") + 
    theme_bw() + 
    facet_wrap(~hour)

ggplot(sp4c_df) + 
    geom_point(aes(x=size, y=max_solar_elevation_difference, color=Latitude), size=5) +
    xlab("Subset size (NxN pixel subset)") + 
    ylab("Max difference in solar elevation between center and four corners") + 
    theme_bw() + 
    facet_wrap(~hour)


##############################################################################
# Do the second thing
##############################################################################

calculate_sp_adj_chunks <- function(tif_filename, date_times, sizes){
    # Load the DEM
    dem <- rast(tif_filename)
    
    # Create lat/long rasters and stack
    lon <- init(dem, 'x')
    lat <- init(dem, 'y')
    s <- c(lon, lat, dem)

    # Calcualte sp in four corners of the subset for each subset size and datetime
    sp_dt = list()
    i = 0
    for(num_pixels in sizes){
        for(d in date_times){
            i = i + 1
            c = round(num_pixels/2)
            solar_position_c = calculate_sp(s[[1:3]][c,c], d[1], d[2]) # solar position in center of top left chunk
            solar_position_below = calculate_sp(s[[1:3]][c+num_pixels,c], d[1], d[2]) # solar position in center of chunk below
            solar_position_nextto = calculate_sp(s[[1:3]][c,c+num_pixels], d[1], d[2]) # sp in center of chunk next to 
            solar_position_kitty = calculate_sp(s[[1:3]][c+num_pixels,c+num_pixels], d[1], d[2]) # sp in center of chunk kittydown      
            dt = rbind(solar_position_below, solar_position_nextto, solar_position_kitty)
            dt$az_diff <- NA
            for(w in 1:nrow(dt)){
                d1 <- abs(dt$solar_azimuth[w] - solar_position_c$solar_azimuth[1])
                d2 <- 360-d1
                dt$az_diff[w] <- min(d1, d2)
            }
            solar_position_c$max_solar_azimuth_difference <- max(dt$az_diff)
            solar_position_c$max_solar_elevation_difference <- max(abs(dt$solar_elevation - solar_position_c$solar_elevation))
            solar_position_c$size = num_pixels
            solar_position_c <- solar_position_c[,!c("solar_azimuth", "solar_elevation")]
            sp_dt[[i]] = solar_position_c
        }
    }
    sp_dt = rbindlist(sp_dt)
    sp_dt$filename <- tif_filename
    return(sp_dt)
}

date_times <- list(c("07", "09"), c("07", "12"), c("07", "15"), c("07", "18"))
sizes <- c(1000, 2000, 3000)
files <- c(
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n00w075_elv.tif",
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n35w085_elv.tif", 
    "C:/Users/kmcquil/Documents/Global_Hillshade/data/raw/merit/n60w075_elv.tif"
    )
out <- "data/outputs/sensitivity_test"

sp_adj_df <- rbindlist(lapply(files, calculate_sp_adj_chunks, date_times=date_times, sizes=sizes))
sp_adj_df$hour <- as.numeric(format(as.POSIXct(sp_adj_df$datetime), format = "%H"))
sp_adj_df$Latitude <- as.numeric(substr(basename(sp_adj_df$filename), 2, 3))
fwrite(sp_adj_df, paste0(out, "/sp_adj_df.csv"))

# Plot max difference between the four corners according to the size 
sp_adj_df$Latitude <- as.factor(sp_adj_df$Latitude)
ggplot(sp_adj_df) + 
    geom_point(aes(x=size, y=max_solar_azimuth_difference, color=Latitude), size=5) +
    xlab("Subset size (NxN pixel subset)") + 
    ylab("Max difference in solar azimuth from other three chunks") + 
    theme_bw() + 
    facet_wrap(~hour)

ggplot(sp_adj_df) + 
    geom_point(aes(x=size, y=max_solar_elevation_difference, color=Latitude), size=5) +
    xlab("Subset size (NxN pixel subset)") + 
    ylab("Max difference in solar elevation from other three chunks") + 
    theme_bw() + 
    facet_wrap(~hour)