###########################################################################
###########################################################################
# Find which tiles + timesteps finished and create a list of the tiles that 
# need to still be run 
###########################################################################
###########################################################################
library(data.table)
# Set wd 
setwd("/projects/swot/kmcquil/Global_Hillshade/")

# Go to the dem tiles folder and get a list of all of the tiles that we 
# wanted to run on 
tiles_list <- list.files("data/dem_tiles", full.names=TRUE)
tiles_dt <- rbindlist(lapply(tiles_list, fread))

# Use the names to create the full filepaths with the shadow type and the date 
outdir <- "data/outputs/shadows_v2"
months <- seq(1, 12)
hours <- seq(0, 23)
dts <- data.table(expand.grid(months, hours))
colnames(dts) <- c("month", "hour")
size <- 2000
all_files <- c()
for(j in 1:nrow(tiles_dt)){
    for(i in 1:nrow(dts)){
        month <- dts$month[i]
        hour <- dts$hour[i]
        tif_filename <- tiles_dt$file[j]
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
        out_fp_raytrace <- paste0(outdir, "/", "raytrace", "/", out_filename)
        out_fp_lambshade <- paste0(outdir, "/", "lambshade", "/", out_filename)
        out_fp_combo <- paste0(outdir, "/", "combo", "/", out_filename)
        all_files <- c(all_files, out_fp_raytrace, out_fp_lambshade, out_fp_combo)
    }
}

# Go to the output folders and get a list of all of the tiles that were finished 
raytrace <- list.files(paste0(outdir, "/", "raytrace"), full.names=TRUE)
lambshade <- list.files(paste0(outdir, "/", "lambshade"), full.names=TRUE)
combo <- list.files(paste0(outdir, "/", "combo"), full.names=TRUE)
processed <- c(raytrace, lambshade, combo)

print(paste0("The number total files to process: ", length(all_files)))
print(paste0("The number of files actually processed so far: ", length(processed)))
print(paste0("fraction completed is: ", length(processed)/length(all_files) ))

# Find the difference and then unique tiles
result <- setdiff(all_files, processed)
unique_tiles <- c()
for(i in result){
    unique_tiles <- c(unique_tiles, substr(basename(i), nchar(basename(i))-14, nchar(basename(i))))
}
unique_tiles <- unique(unique_tiles)

fps <- c()
for(i in unique_tiles){
    fps <- c(fps, 
        paste0(
        dirname(tiles_dt$file[1]), 
        "/", 
        i
        )
    )  
}

# Create new files
n <- 30
r <- ceiling(length(fps) / n)
b <- seq(1, length(fps), by = r)
b <- c(b, length(fps) + 1)


# Make the output directory if it doesn't already exist 
outdir <- "data/dem_tiles_cleanup"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
unlink(paste0(outdir, "/*"), recursive = FALSE)

for (i in seq_along(b)[-length(b)]) {
  sub <- fps[b[i]:(b[i+1] - 1)]
  df <- data.frame(file = unlist(sub))
  j <- i
  print(j)
  write.csv(df, file = file.path(outdir, paste0("tiles", j, ".csv")), 
         row.names = FALSE, quote = TRUE)
}