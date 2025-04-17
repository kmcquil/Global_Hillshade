import itertools
import pandas as pd
import math

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
fps = []
for i in range(0, len(combinations)):
    fps.append("data/raw/merit_retile/" + combinations[i][0] + combinations[i][1] + "_elv.tif")

fps = pd.DataFrame(fps)
fps.columns = ["file"]

r = math.ceil(len(fps)/10)
b = list(range(0, len(fps) ,r))
b.append(len(fps))

for i in range(0, 10):
    sub = fps.iloc[b[i]:b[i+1]]
    j = i + 1
    sub.to_csv("data/tiles_conus/tiles"+str(j) + ".csv", index=False)