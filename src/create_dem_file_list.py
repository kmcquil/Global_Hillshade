import pandas as pd
import itertools
import math
import glob 
import os

# List all of the retiled DEM tifs and make into a pd df
fps = glob.glob("data/raw/merit_retile/*.tif")
fps = pd.DataFrame(fps)
fps.columns = ["file"]

# Don't include the tiles already processed from conus
lats = list(range(25, 50, 5))
lats = ["n" + str(i).zfill(2) for i in lats]
longs = list(range(70, 130, 5))
longs = ["w" + str(i).zfill(3) for i in longs]
combinations = list(itertools.product(lats, longs))
all_fps = []
for i in range(0, len(combinations)):
    all_fps.append("data/raw/merit_retile/" + combinations[i][0] + combinations[i][1] + "_elv.tif")
fps = fps[~fps['file'].isin(all_fps)]

# Check if they exist
fps = fps[fps['file'].apply(os.path.exists)]
print("Number of files: " + str(len(fps)))

# Make sequence of fps to spread across n files
n = 30
r = math.ceil(len(fps)/n)
b = list(range(0, len(fps) ,r))
b.append(len(fps))

# Make the output directory if it doesn't already exist 
outdir = "data/dem_tiles"
if not os.path.exists(outdir):
    os.makedirs(outdir)

for i in range(0, len(b)-1):
    sub = fps.iloc[b[i]:b[i+1]]
    j = i + 1
    sub.to_csv(outdir+"/tiles"+str(j) + ".csv", index=False)