# Test how big files will be based on saving in different ways 
import rasterio
import numpy as np

fp = r"C:\Users\kmcquil\Documents\Global_Hillshade\data\outputs\scratch\shadow_2000_7_15_n40w100_elv.tif"

###################################################################
# Just changing to int 
###################################################################
# Step 1: Open the source .tif
with rasterio.open(fp) as src:
    float_data = src.read(1)  # Read the first band (or loop for multiband)
    profile = src.profile     # Copy profile for later
    transform = src.transform

# Convert no data value to integer 
# The new no data value will be set to 9999
# So here make it 9999*1000 since we divide by that in the next step
float_data[np.isnan(float_data)] = 9999000

# Step 2: Scale floats by 1000 and convert to uint16
scaled_data = (float_data * 1000).round().astype('uint16')

# Step 3: Update profile to match new dtype
profile.update(
    dtype='uint16',
    nodata=9999,
    #compress='LZW',
    #tiled=True,
    #blockxsize=256,
    #blockysize=256
)

# Step 4: Save the new raster
out = r"C:\Users\kmcquil\Documents\Global_Hillshade\data\outputs\scratch\just_int.tif"
with rasterio.open(out, 'w', **profile) as dst:
    dst.write(scaled_data, 1)  # Write to band 1

###################################################################
# Int and compression
###################################################################
# Step 1: Open the source .tif
with rasterio.open(fp) as src:
    float_data = src.read(1)  # Read the first band (or loop for multiband)
    profile = src.profile     # Copy profile for later
    transform = src.transform

# Convert no data value to integer 
# The new no data value will be set to 9999
# So here make it 9999*1000 since we divide by that in the next step
float_data[np.isnan(float_data)] = 9999000

# Step 2: Scale floats by 1000 and convert to uint16
scaled_data = (float_data * 1000).round().astype('uint16')

# Step 3: Update profile to match new dtype
profile.update(
    dtype='uint16',
    nodata=9999,
    compress='LZW',
    #tiled=True,
    #blockxsize=256,
    #blockysize=256
)

# Step 4: Save the new raster
out = r"C:\Users\kmcquil\Documents\Global_Hillshade\data\outputs\scratch\int_and_compress.tif"
with rasterio.open(out, 'w', **profile) as dst:
    dst.write(scaled_data, 1)  # Write to band 1


###################################################################
# Just compression
###################################################################
# Step 1: Open the source .tif
with rasterio.open(fp) as src:
    float_data = src.read(1)  # Read the first band (or loop for multiband)
    profile = src.profile     # Copy profile for later
    transform = src.transform

# Step 3: Update profile to match new dtype
profile.update(
    compress='LZW',
    #tiled=True,
    #blockxsize=256,
    #blockysize=256
)

# Step 4: Save the new raster
out = r"C:\Users\kmcquil\Documents\Global_Hillshade\data\outputs\scratch\just_compress.tif"
with rasterio.open(out, 'w', **profile) as dst:
    dst.write(float_data, 1)  # Write to band 1


###################################################################
# try uint8
###################################################################
# Step 1: Open the source .tif
with rasterio.open(fp) as src:
    float_data = src.read(1)  # Read the first band (or loop for multiband)
    profile = src.profile     # Copy profile for later
    transform = src.transform

# Convert no data value to integer 
# The new no data value will be set to 255
# 255 = x * 100
# So here make it 255*100 since we divide by that in the next step
float_data[np.isnan(float_data)] = 25500

# Step 2: Scale floats by 1000 and convert to uint16
scaled_data = (float_data * 100).round().astype('uint8')

# Step 3: Update profile to match new dtype
profile.update(
    dtype='uint8',
    nodata=255,
    compress='LZW',
    #tiled=True,
    #blockxsize=256,
    #blockysize=256
)

# Step 4: Save the new raster
out = r"C:\Users\kmcquil\Documents\Global_Hillshade\data\outputs\scratch\int_and_compress_uint8.tif"
with rasterio.open(out, 'w', **profile) as dst:
    dst.write(scaled_data, 1)  # Write to band 1

