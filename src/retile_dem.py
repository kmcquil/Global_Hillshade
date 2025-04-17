##############################################################################
##############################################################################
# Retile the DEM to account for edge effects and long shadows
##############################################################################
##############################################################################

# Concept: For each tile, we calculate the longest possible shadow and add that 
# length to the bounds of the raster to create a bigger raster. By doing this for 
# each tile, we create a mosaic of overlapping tiles which should account for 
# long shadows and edge effects 

##############################################################################
# Setup
##############################################################################
import tarfile 
from osgeo import gdal
import glob
import math
import rioxarray
import geopandas as gpd
from shapely.geometry import box
import pandas as pd
import os
import rasterio
from rasterio.mask import mask

##############################################################################
# Functions
##############################################################################

def extract_tar(tar_filepath, dest_dir):
    """
    Extracts a tar archive to a specified directory.
    Args:
        tar_filepath (str): The path to the tar archive file.
        dest_dir (str): The directory to extract the contents to.
    Return (null)
    """
    #with tarfile.open(tar_filepath, 'r') as tar:
    #    tar.extractall(dest_dir)
    with tarfile.open(tar_filepath, 'r') as tar_file:
        for member in tar_file.getmembers():
            if member.isfile():
                member.name = os.path.basename(member.name)
                tar_file.extract(member, path=dest_dir)

def shadow_length(elevation, angle):
    """
    Calculates the maximum shadow length for a given dem
    Args:
        elevation (float): elevation (meters)
        angle (float): angle to calculate max shadow length at
    Return (float) maximum shadow length in meters
    """
    dist_meters = elevation/math.tan((angle*math.pi)/180)
    return dist_meters

def convert_wgs_to_utm(lon, lat):
    """
    Identify the correct utm zone for epsg4326 coordinates
    Args:
        lon (float): longitude in 4326
        lat (float): latitude in 4326
    Return (integer) epsg code for the utm zone
    """
    utm_band = str(math.floor((lon+180)/6)+1)
    utm_band = utm_band.zfill(2)
    if lat>0:
        epsg_code = '326'+utm_band
        return epsg_code
    else:
        epsg_code = '327'+utm_band
        return epsg_code

def calculate_new_tile_bounds(tile):
    """
    Calculate new tile boundaries based on the max shadow length
    Args:
        tile (string): filepath to the tile 
    Return (pd df) includes info on the old and new tile bounds
    """
    with rioxarray.open_rasterio(tile) as ds:
        max_elevation = float(ds.max())
        max_shadow_length = shadow_length(max_elevation, 5)

        # Convert raster bounds to geopandas
        xmin_og, ymin_og, xmax_og, ymax_og = ds.rio.bounds()

    bbox = box(xmin_og, ymin_og, xmax_og, ymax_og)
    geo_df = gpd.GeoDataFrame({'geometry': [bbox]}, crs=ds.rio.crs)

    # Identify the correct utm crs and reproject 
    utm_crs = convert_wgs_to_utm((xmin_og+xmax_og)/2, (ymin_og+ymax_og)/2)
    geo_df = geo_df.to_crs(utm_crs)

    # Get the new bounds and add the shadow length
    xmin, ymin, xmax, ymax = geo_df.total_bounds
    xmin = xmin-max_shadow_length
    xmax = xmax+max_shadow_length
    ymin = ymin-max_shadow_length
    ymax = ymax+max_shadow_length

    # Put into a gdf again to convert back to wgs84 and get the final bounds
    bbox = box(xmin, ymin, xmax, ymax)
    geo_df = gpd.GeoDataFrame({'geometry': [bbox]}, crs=geo_df.crs)
    geo_df = geo_df.to_crs(4326)
    xmin, ymin, xmax, ymax = geo_df.total_bounds
    
    # Deal with edge cases where the tile is trying to wrap around the world 
    # bc it goes past -180 or 180 to the other side
    if xmax_og > 175: 
        xmin = xmax.copy() 
        xmax = 180
    if xmax_og < -175: 
        xmax = xmin.copy() 
        xmin = -180

    # Save these to a df 
    df = pd.DataFrame({
        'tile':[os.path.basename(tile)],
        'max_elevation':[max_elevation],
        'max_shadow_length':[max_shadow_length],
        'xmin_og':[xmin_og],
        'xmax_og':[xmax_og],
        'ymin_og':[ymin_og],
        'ymax_og':[ymax_og],
        'xmin_rt':[xmin],
        'xmax_rt':[xmax],
        'ymin_rt':[ymin],
        'ymax_rt':[ymax]
    })
    return df

def clip_vrt_to_coords(vrt_path, output_path, min_x, min_y, max_x, max_y):
    """
    Clips a VRT file to the specified coordinates and saves the result.
    Args:
        vrt_path (str): Path to the input VRT file.
        output_path (str): Path to save the clipped raster.
        min_x (float): Minimum X coordinate.
        min_y (float): Minimum Y coordinate.
        max_x (float): Maximum X coordinate.
        max_y (float): Maximum Y coordinate.
    Return (null) write out the retiled tile
    """
    with rasterio.open(vrt_path) as src:
        # Create a bounding box geometry from coordinates
        bbox = box(min_x, min_y, max_x, max_y)
        
        # Transform the bounding box to the VRT's coordinate reference system (CRS)
        # This is just to get it into the right format for rasterio mask. 
        #geom = rasterio.warp.transform_geom('EPSG:4326', src.crs, bbox.__geo_interface__)
        geom = rasterio.warp.transform_geom(src.crs, 
                                            src.crs, 
                                            bbox.__geo_interface__)
        # Clip the raster using the transformed geometry
        out_image, out_transform = mask(src, [geom], crop=True)
        out_meta = src.meta.copy()

        # Update metadata for the clipped raster
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        # Save the clipped raster
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image)


##############################################################################
# Unzip each .tar file
##############################################################################
# Set the path to extract files
extraction_path = "data/raw/merit"
# List the tar files
input_files = glob.glob("data/raw/merit/*.tar")
# Extract
[extract_tar(file, extraction_path) for file in input_files]

##############################################################################
# Create a virtual raster of all of the tiles
##############################################################################
input_files = glob.glob("data/raw/merit/*.tif", recursive=True)
output_vrt = "data/raw/merit_retile/elv.vrt"
vrt = gdal.BuildVRT(output_vrt, input_files)
del vrt

##############################################################################
# Calculate new tile bounds based on the potential maximum shadow length
##############################################################################
tile_bounds = pd.concat([calculate_new_tile_bounds(tile) for tile in input_files])
tile_bounds.to_csv("data/outputs/tile_bounds.csv")

##############################################################################
# Create new tiles from the VRT
##############################################################################
vrt_file = "data/raw/merit_retile/elv.vrt"
for i in range(0, tile_bounds.shape[0]):
    output_file = "data/raw/merit_retile/" + tile_bounds['tile'].iloc[i]
    min_x_coord = tile_bounds['xmin_rt'].iloc[i]
    min_y_coord = tile_bounds['ymin_rt'].iloc[i]
    max_x_coord = tile_bounds['xmax_rt'].iloc[i]
    max_y_coord = tile_bounds['ymax_rt'].iloc[i]
    clip_vrt_to_coords(vrt_file, output_file, min_x_coord, min_y_coord, max_x_coord, max_y_coord)

##############################################################################
# Delete the un-tarred tif files 
##############################################################################
input_files = glob.glob("data/raw/merit/*.tif")
[os.remove(file) for file in input_files]

##############################################################################
# For visualization purposes, make a geopandas df of the original and 
# retiled tiles 
##############################################################################
tile_gdf = []
for i in range(0, tile_bounds.shape[0]):
    bbox = box(tile_bounds['xmin_og'].iloc[i], 
               tile_bounds['ymin_og'].iloc[i], 
               tile_bounds['xmax_og'].iloc[i],
               tile_bounds['ymax_og'].iloc[i])
    gdf = gpd.GeoDataFrame({
        'geometry':[bbox],
        'version':['Original'],
        'tile':[tile_bounds['tile'].iloc[i]]
        }, crs=4326)
    tile_gdf.append(gdf)

    bbox = box(tile_bounds['xmin_rt'].iloc[i], 
               tile_bounds['ymin_rt'].iloc[i], 
               tile_bounds['xmax_rt'].iloc[i],
               tile_bounds['ymax_rt'].iloc[i])
    gdf = gpd.GeoDataFrame({
        'geometry':[bbox],
        'version':['Retiled'],
        'tile':[tile_bounds['tile'].iloc[i]]
        }, crs=4326)
    tile_gdf.append(gdf)
tile_gdf =pd.concat(tile_gdf)
tile_gdf.to_file("data/outputs/tile_bounds.shp")

