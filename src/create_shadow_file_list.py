import pandas as pd
import math
import glob 
import os

# List all of the retiled DEM tifs and make into a pd df
fps = glob.glob("data/outputs/shadows/*.tif")
fps = pd.DataFrame(fps)
fps.columns = ["file"]

# Make sequence of fps to spread across n files
n = 30
r = math.ceil(len(fps)/n)
b = list(range(0, len(fps) ,r))
b.append(len(fps))

# Make the output directory if it doesn't already exist 
outdir = "data/shadow_tiles"
if not os.path.exists(outdir):
    os.makedirs(outdir)

for i in range(0, n):
    sub = fps.iloc[b[i]:b[i+1]]
    j = i + 1
    sub.to_csv(outdir+"/tiles"+str(j) + ".csv", index=False)