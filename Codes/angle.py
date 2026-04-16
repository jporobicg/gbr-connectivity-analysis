## File: angle.py
## date : 2022-11-08
## creator: Javier Porobic
## email: javier.porobicgarate@csiro.au
## Description : calculates the angle between the centroids of 2 reefs and creates a csv file with the angle value
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~ ##
## ~     Libraries  ~ ##
## ~~~~~~~~~~~~~~~~~~ ##
import geopandas as gpd
import numpy as np
import math



## ~~~~~~~~~~~~~~~~~~ ##
## ~     Functions  ~ ##
## ~~~~~~~~~~~~~~~~~~ ##
def veclength(vector):
    """This function calculates the size (length) of a vector. in only needs the vector of size 2
    """
    value = math.sqrt(math.pow(vector[0], 2) + math.pow(vector[1], 2))
    return value

def angle(a,b):
    """ this function calculates the angle between two reefs
    I should update this function using the bearing formula
    """
    dp = np.dot(a, b) ## Dot product of the two vectors
    la = veclength(a)
    lb = veclength(b)
    costheta =dp / (la * lb)
    rads = math.acos(costheta)
    angle = 180 * rads / math.pi
    return (angle)

def haversine(coord1, coord2):
    # Earth's radius in km
    radius = 6371
    # Convert coordinates to radians
    lat1, lon1 = np.radians(coord1)
    lat2, lon2 = np.radians(coord2)
    # Calculate differences between coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    # Calculate Haversine formula
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    # Calculate distance in km
    distance = radius * c
    return distance




## ~~~~~~~~~~ ##
## ~   main ~ ##
## ~~~~~~~~~~ ##
shapefile = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/reefs/AIMS_shapefile/Simplified_LTMS_AIMS.shp'
#shapefile = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/Kitchen_model/test_shapefile.shp'
data_shape = gpd.read_file(shapefile)
num_sites = data_shape.shape[0]
reef_centroids = []
## getting the centroid's location
for site in range(0, num_sites):
    ## release reef
    value_index = list(data_shape.loc[data_shape['FID'] == site].index)
    value_index = int("".join(map(str, value_index)))
    polygon = data_shape['geometry'][value_index]
    reef_centroids.append(polygon.centroid)


## creating the output file
direction_outFile = open('/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/outputs/GBR_zones_direction.csv', 'w')
angle_outFile = open('/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/outputs/GBR_zones_angles.csv', 'w')
distance_outFile = open('/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/outputs/GBR_zones_distance.csv', 'w')

## adding the headers
angle_outFile.write("ReleaseSite\TargetSite")
direction_outFile.write("ReleaseSite\TargetSite")
distance_outFile.write("ReleaseSite\TargetSite")
for site in range(0,num_sites):
    value_index = list(data_shape.loc[data_shape['FID'] == site].index)
    value_index = str("".join(map(str, value_index)))
    angle_outFile.write(',' + value_index)
    direction_outFile.write(',' + value_index)
    distance_outFile.write(','  + value_index)


angle_outFile.write("\n")
direction_outFile.write("\n")
distance_outFile.write("\n")
## vecot along the Y axis angle 0
a = [0, 1]
for release_site in range(0, num_sites):
    value_index = list(data_shape.loc[data_shape['FID'] == release_site].index)
    value_index = str("".join(map(str, value_index)))
    angle_outFile.write(str(value_index))
    direction_outFile.write(str(value_index))
    distance_outFile.write(str(value_index))
    for target_site in range(0, num_sites):
        """This will do the calculations for each vector.
        Normalized distance between the two reefs.
        """
        reef_angle = 0
        direction = 0
        distance = 0
        if release_site != target_site:
            coordinates_sink = np.array([reef_centroids[target_site].coords[0][0], reef_centroids[target_site].coords[0][1]])
            coordinates_source =  np.array([reef_centroids[release_site].coords[0][0], reef_centroids[release_site].coords[0][1]])
            b =  [reef_centroids[target_site].coords[0][0] - reef_centroids[release_site].coords[0][0], reef_centroids[target_site].coords[0][1] - reef_centroids[release_site].coords[0][1] ]
            reef_angle = angle(a,b)
            if b[0] < 0:
                reef_angle = 360 - reef_angle
            ##Adding the rotation reef_angle from the hydro model to the shape file
            rot_reef_angle = reef_angle + 22.5
            ## asign the sector (every 10deg so 36 sector) to each reef reef reef_angle
            direction = math.floor(rot_reef_angle/ 10) % 36
            distance = haversine(coordinates_source, coordinates_sink)
        angle_outFile.write(',' + str(reef_angle))
        direction_outFile.write(',' + str(direction))
        distance_outFile.write(',' + str(distance))
    angle_outFile.write("\n")
    direction_outFile.write("\n")
    distance_outFile.write("\n")

angle_outFile.close()
direction_outFile.close()
distance_outFile.close()
