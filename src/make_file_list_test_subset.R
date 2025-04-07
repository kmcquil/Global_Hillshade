##############################################################################
# Create lists of tiles to use job array for parallel shadow calculations
##############################################################################

# Load libraries
library(data.table)

# Create a folder to the tile csvs
dir.create("data/tiles", showWarnings = FALSE)

# List the files 
files <- data.table(file=c("data/raw/merit_retile/n35w115_elv.tif", "data/raw/merit_retile/n35w110_elv.tif", "data/raw/merit_retile/n40w115_elv.tif", "data/raw/merit_retile/n40w110_elv.tif"))

# Choose the number of chunks
n <- 4
#breaks <- round(seq(from=1, to=nrow(files), length.out=(n+1)))
# Save as csvs
#for(i in 1:n){
    #dt <- files[breaks[i]:(breaks[i+1]-1),]
    #fwrite(dt[,"file"], paste0("data/tiles/tiles", i, ".csv"))
#}
breaks <- round(seq(from=1, to=nrow(files)))
for(i in 1:n){
    dt <- files[breaks[i]:breaks[i],]
    fwrite(dt[,"file"], paste0("data/tiles/tiles", i, ".csv"))
}
