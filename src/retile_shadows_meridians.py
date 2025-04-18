##############################################################################
# Mosaic the newly mapped shadows and convert back to original tile sizes
##############################################################################

##############################################################################
# Load modules
##############################################################################
import glob
import pandas as pd
import os
import rasterio
from rasterio.merge import merge
import numpy as np

##############################################################################
# Retile 
##############################################################################

# Open the df with original tile bounds 
tile_bounds = pd.read_csv("data/outputs/tile_bounds.csv")

def custom_merge_avg(old_data, new_data, old_nodata, new_nodata, index=None, roff=None, coff=None):
    old_data[:] = np.nanmean( np.array([ old_data, new_data ]), axis=0 )

def mosaic_and_retile(file, outdir, merge_method):
    # Check if the file exists, if not return 
    if not os.path.exists(file):
        print("doesn't exist")
        return

    # Get the original bounds for this tile from the tile bounds df 
    dem_tile_name = os.path.basename(file)[-15:]
    tile_bounds_subset = tile_bounds[tile_bounds["tile"]==dem_tile_name]
    # (left, bottom, right, top)
    bbox = (tile_bounds_subset['xmin_og'].iloc[0],
               tile_bounds_subset['ymin_og'].iloc[0],
               tile_bounds_subset['xmax_og'].iloc[0],
               tile_bounds_subset['ymax_og'].iloc[0])

    # Extract the lat and long 
    file_split = os.path.basename(file).split("_")
    lat_ = file_split[4][0:3]
    long_ = file_split[4][3:7]

    # Identify the tiles that should be up/down and left/right -- account for edge cases on N/S and E/W
    # This tiles are in 5 degree increments
    # For latitude, N starts at 0 and S starts at 5
    # For longitude, E starts at 0 and W starts at 
    
    # Start with latitude
    lat_dir = lat_[0]
    lat_num = lat_[1:]
    lat_num_south = str(int(lat_num)-5).zfill(2)
    lat_dir_south = lat_dir
    lat_num_north = str(int(lat_num)+5).zfill(2)
    lat_dir_north= lat_dir     

    # Edge case at N-S divide on the N side
    if ((lat_num == '00') & (lat_dir == 'n')):
        lat_num_south = '05'
        lat_dir_south = 's'
        lat_num_north = '05'
        lat_dir_north = 'n'
    # Edge case at N-S divide on the S side 
    if ((lat_num == '05') & (lat_dir == 's')):
        lat_num_south = '10'
        lat_dir_south = 's'
        lat_num_north = '00'
        lat_dir_north =  'n'          
    
    # Now longitude
    long_dir = long_[0]
    long_num = long_[1:]
    long_num_east = str(int(long_num)-5).zfill(3)
    long_dir_east = long_dir
    long_num_west = str(int(long_num)+5).zfill(3)
    long_dir_west = long_dir   

    # edge case of E/W divide on E side
    if ((long_num == '000') & (long_dir == 'e')):
        long_num_east = '005'
        long_dir_east = 'e'
        long_num_west = '005'
        long_dir_west = 'w'
    # edge case of E/W divide on E side 
    if ((long_num == '005') & (long_dir == 'w')):
        long_num_east = '000'
        long_dir_east = 'e'
        long_num_west = '010'
        long_dir_west = 'w'

    # Combine all of this into filenames 
    # up/down file is the same longitude but new latitudes 
    file_up = file_split.copy()
    file_up[4] = lat_dir_north + lat_num_north + long_dir + long_num
    file_up = "_".join(file_up)
    file_down = file_split.copy()
    file_down[4] = lat_dir_south + lat_num_south + long_dir + long_num
    file_down = "_".join(file_down)
    file_east = file_split.copy()
    file_east[4] = lat_dir + lat_num + long_dir_east + long_num_east
    file_east = "_".join(file_east)
    file_west = file_split.copy()
    file_west[4] = lat_dir + lat_num + long_dir_west + long_num_west
    file_west = "_".join(file_west)
    adj_files_ = [file_up, file_down, file_east, file_west]
    adj_files_ = [os.path.dirname(file) + "/" + i for i in adj_files_]
    adj_files_.append(file)

    # Check they actually exist
    adj_files = []
    for i in adj_files_:
        if os.path.exists(i):
            adj_files.append(i)
    print(len(adj_files))

    # Mosaic the files together, using the bounds of the original tile 
    src_files_to_mosaic = [rasterio.open(s) for s in adj_files]
    mosaic, out_trans = merge(src_files_to_mosaic, bounds=bbox, method=merge_method, nodata=float('nan')) 
    
    # Update the metadata with new dimensions, transform, and CRS
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": src_files_to_mosaic[0].crs
    })

    # Write the mosaic raster to a new file
    output_path = outdir + "/" + os.path.basename(file)
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    # Close the opened raster datasets
    [src.close() for src in src_files_to_mosaic]
    print("finished")

# List all of the shadow-mapped files
shadow_files = glob.glob("data/outputs/test_meridians/*.tif")
outdir = "data/outputs/test_meridians_retile_max"
if not os.path.exists(outdir):
    os.makedirs(outdir)
[mosaic_and_retile(file, outdir, 'max') for file in shadow_files]

print("finished the retiling")

def big_mosaic(file_list, output_path, merge_method):
    src_files_to_mosaic = [rasterio.open(s) for s in file_list]
    mosaic, out_trans = merge(src_files_to_mosaic, method=merge_method) 
    
    # Update the metadata with new dimensions, transform, and CRS
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": src_files_to_mosaic[0].crs
    })

    # Write the mosaic raster to a new file
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    # Close the opened raster datasets
    [src.close() for src in src_files_to_mosaic]
    print("finished")


# Create list of unique date times
shadow_files = glob.glob(outdir + "/*.tif")
def get_dt(file):
    z = os.path.basename(file).split("_")[0:4]
    z = "_".join(z)
    return z
dts = list(set([get_dt(file) for file in shadow_files]))
def big_mosaic_dt(dt):
    print(dt)
    files = glob.glob(outdir + "/" + dt + "*.tif")
    big_mosaic(files, outdir + "/all_together_" + dt + ".tif", "max")
[big_mosaic_dt(dt) for dt in dts]