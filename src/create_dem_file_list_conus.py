import itertools
import pandas as pd
import math
import os

# Range for conus
#25N - 45N 
#70W - 125W 

# Create the combinations 
lats = list(range(25, 50, 5))
lats = ["n" + str(i).zfill(2) for i in lats]

longs = list(range(70, 130, 5))
longs = ["w" + str(i).zfill(3) for i in longs]

combinations = list(itertools.product(lats, longs))

# Create the file paths
all_fps = []
for i in range(0, len(combinations)):
    all_fps.append("data/raw/merit_retile/" + combinations[i][0] + combinations[i][1] + "_elv.tif")

# Check if they exist
fps = [i for i in all_fps if os.path.exists(i)]
fps = pd.DataFrame(fps)
fps.columns = ["file"]
print("Number of files: " + str(len(fps)))

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
    print(j)
    sub.to_csv(outdir+"/tiles"+str(j) + ".csv", index=False)