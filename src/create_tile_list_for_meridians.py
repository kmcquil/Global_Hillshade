import itertools
import pandas as pd
import math

# Range for meridians
#5S - 5N
#5E - 10E

# Create the combinations 
lats = ["s05", "n00"]
longs = ["e005", "e010"]
combinations = list(itertools.product(lats, longs))

# Create the file paths
fps = []
for i in range(0, len(combinations)):
    fps.append("data/raw/merit_retile/" + combinations[i][0] + combinations[i][1] + "_elv.tif")

fps = pd.DataFrame(fps)
fps.columns = ["file"]

r = math.ceil(len(fps)/4)
b = list(range(0, len(fps) ,r))
b.append(len(fps))

for i in range(0, 4):
    sub = fps.iloc[b[i]:b[i+1]]
    j = i + 1
    sub.to_csv("data/tiles_meridian/tiles"+str(j) + ".csv", index=False)