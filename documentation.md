# Map Shadows

## Data sources
- Merit Hydro DEM: https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/

## Rayshader 
- I updated the source code in the Rayshader package. 
    - https://github.com/kmcquil/rayshader_km
    - Return separate outputs for lamb shading, ray shading, and the combined product
    - Stop rescaling shadow intensity from 0-1 and instead enforce a 0 min and 1 max. The issue with the rescaling was that it made adjacent tiles inconsistent.
- Ray tracing is performed for each pixel by supplying the solar angle for the pixel location based on the date-time. Therefore, it is most accurate by calcualting the solar angle for each pixel and then perfomring ray tracing. But that is really slow. To speed up calculations, calculate the solar position at the center of 2000 x 2000 pixel blocks and use that for all 2000x2000 pixels. Each tif has ~6000x6000 pixels, therefore solar position is calculated and applied for ~9 subsets. This isn't perfect but attempts to balance precision and speed.

## Compute environment
- The rayshader package was difficult to install directly onto HPC so instead I created a docker/sif container to reliably use the package on HPC. See docker/notes.md for details on how this worked. 

## Scripts
Each script has a corresponding .sh file that is named the same with _submit.sh added to the end.

- data_download.py 
    - Script downloads the merit hydro dem data
- retile_dem
    - Script retiles the merit hydro dem to account for shadows across tiles. For each tile, calculate the maximum possible shadow length based on topograhpy and add that distance to all sides of tile.
- create_dem_file_list.py / create_dem_file_list_conus.py
    - Script creates csvs with lists of retiled merit dems to parallelize shadow mapping
- map_shadows.R
    - Script maps shadows using map_shadows_submit.sh. Due to HPC time/memory constraints, all of the shadow mapping can't be done in one submission.  
    - After this script is run for the first time, run tracking.R. Tracking.R checks how many of the files successfully completed and then creates new csvs with lists of retiled merit dems that shadows still need to be calcualted for. 
    - map_shadows_cleaing_submit.sh starts shadow calculations on the list of remianing tiles from the prior step. Continue this step and last step until all tiles are completed. 
- create_shadow_file_list.py 
    - Script creates csvs with lists of shadow tifs that need to be retiled in parallel
- retile_shadows.py
    - retile the shadow tifs to the original size/shape of the merit dem. In areas where the tiles overlap, take the max.
- sensitivity.R
    - This script was run before the global push to understand the sensitivity of shadow intensity and solar position with different size chunks of pixels within a tif to guide our selection of precision/speed.

## Data
- raw
    - merit: Merit Hydro DEM tiles
    - merti_retile: Retiled Merit Hydro DEM tiles. These are used in shadow mapping.
    - sensitivity_examples: Merit Hydro DEM tiles that were used for sensitivity analyses of the solar position and shadow intensity.
- dem_tiles: Contains csvs with lists of DEM tiles. Used to submit jobs in parallel.
- outputs
    - scratch: old stuff i was scared to delete before.
    - sensitivity_test: results from sensitivity test
    - shadows: Outputs of shadow mapping broken into three separate products
        - raytrace: Shadows from raytracing
        - lambshade: Shadows from lambert shading
        - combo: Combined product of lambert shading and ray tracing
    - shadows_retile: Shadow mapping outputs retiled to the original merit hydro dem tiles. This is the final product.

## HPC
- This project can be found at "/projects/swot/kmcquil/Global_Hillshade/"
- I added Yohtaro to the project, but if there are permissions issues that come up we can trouble shoot.

## Outputs
- Rayshader produces maps of shadow intensity from 0-1. To save the data more efficiently as Unsigned 8 bit integers, NA values were converted to 255 and shadow intensity was scaled by 100 to range from 0 - 100.