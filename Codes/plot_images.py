import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

# Load the shapefile
gdf = gpd.read_file('../datasets/shapefiles/gbr1_coral_1m_merged.shp')
gdf = gdf.to_crs("EPSG:3857")  # Convert to Web Mercator
# Calculate the centroids once and store them in a dictionary
centroids = gdf.geometry.centroid
centroids_dict = {fid: centroid for fid, centroid in zip(gdf['FID'], centroids)}


# Load the connectivity matrix
#df = pd.read_csv('../datasets/outputs/2015-10-29_T_Connectivity_max.csv', header=None)
df = pd.read_csv('../datasets/outputs/2015-10-29_Connectivity_max.csv', header=None)
# Create a dictionary
# Select the first 100 rows of the DataFrame
subset = df
#  Convert the subset to a dictionary
subset_dict = subset.to_dict('index')

# Create a plot
fig, ax = plt.subplots(figsize=(10, 10))  # Adjust the size as needed

# Plot the shapefile
gdf.plot(ax=ax)

# Set the x and y limits to the top left corner
#ax.set_xlim(left=142, right=145)  # replace with your coordinates
#ax.set_ylim(bottom=-12, top=-10)  # replace with your coordinates
# Iterate over the dictionary
# Remove items with a value of 0
subset_dict = {source: {destination: value for destination, value in connections.items() if value != 0} 
               for source, connections in subset_dict.items()}

# for source, connections in subset_dict.items():
#     for destination, value in connections.items():
#         print(source, destination, value)
# # Iterate over the dictionary
for source, connections in subset_dict.items():
    for destination, value in connections.items():
        if value > 0:  # Assuming you want to plot only if there's a connection
            source_coords = centroids_dict[source]
            destination_coords = centroids_dict[destination]
            plt.plot([source_coords.x, destination_coords.x], [source_coords.y, destination_coords.y], color='black', alpha=value*5)

plt.show()
## create a png file
fig.savefig('output.png', dpi=300)
