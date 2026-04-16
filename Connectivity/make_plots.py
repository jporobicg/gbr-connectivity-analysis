# float connectivity(time, source, sink, treatment, sample)

import geopandas as gpd
import xarray as xr

from plot_images_functions import plot_probability_of_connection, heatmap_connectivity_by_year

gdf = gpd.read_file('/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/shapefiles/AIMS_shapefile/Simplified_LTMS_AIMS.shp')
nc_dataset = xr.open_dataset('/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/connectivity_matrices/connectivity_merulinidae_single.nc')
output_folder = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/figures/'
species = 'merulinidae'
plot_probability_of_connection(nc_dataset, gdf, output_folder + f'probability_of_connection_{species}.png')


heatmap_connectivity_by_year(nc_dataset, output_folder + f'heatmap_connectivity_by_year_{species}.png')
nc_dataset.close()