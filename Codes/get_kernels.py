import geopandas as gpd
import numpy as np
import math
import glob
import xarray as xr


# Loading Functions
# Load input files for connectivity, angle, and shape information
def get_centroids(shapefile):
    data_shape = gpd.read_file(shapefile)
    num_sites = data_shape.shape[0]
    reef_centroids = []
    reef_order = []
    # getting the centroid's location
    for site in range(0, num_sites):
        # release reef
        value_index = list(data_shape.loc[data_shape['FID'] == site].index)
        value_index = int("".join(map(str, value_index)))
        reef_order.append(value_index)
        polygon = data_shape['geometry'][value_index]
        reef_centroids.append(polygon.centroid)
    return reef_centroids, reef_order


def reef_ordering(distance_matrix, angle_matrix, reef_order):
    """""This function will order the outputs using the order that is coming from the shapefile"""
    idx = np.empty_like(reef_order)
    idx[reef_order] = np.arange(len(reef_order))
    distance_matrix_ordered = distance_matrix[:, idx][idx]
    angle_matrix_ordered = angle_matrix[:, idx][idx]
    return distance_matrix_ordered, angle_matrix_ordered


def select_source_locations(connectivity_matrix, angle_matrix, distance_matrix, sink_area):
    # Select source locations that contribute to recruitment in the sink area
    # based on the connectivity values
    source_connection = connectivity_matrix[:, sink_area]
    source_angles = angle_matrix[:, sink_area]
    source_distance = distance_matrix[:, sink_area]
    return source_connection, source_angles, source_distance


def select_sink_locations(connectivity_matrix, angle_matrix, distance_matrix, source_area):
    # Select sink locations that contribute to recruitment in the sink area
    # based on the connectivity values
    sink_connection = connectivity_matrix[source_area, :]
    sink_angles = angle_matrix[source_area, :]
    sink_distance = distance_matrix[source_area, :]
    return sink_connection, sink_angles, sink_distance


def remove_outliers(connectivity, angles, centroids, distance, percent):
    # Remove the percent of locations with smallest and largest recruitment values
    num_outliers = int(len(np.where(connectivity)[0] > 0) * percent / 100)
    sorted_indices = np.argsort(connectivity)
    if (num_outliers == 0):
        filtered_indices = sorted_indices
    else:
        filtered_indices = sorted_indices[num_outliers:-num_outliers]
    filtered_connectivity = [connectivity[i] for i in filtered_indices]
    filtered_distance = [distance[i] for i in filtered_indices]
    filtered_angles = [angles[i] for i in filtered_indices]
    filtered_centroids = [centroids[i] for i in filtered_indices]
    return filtered_connectivity, filtered_angles, filtered_centroids, filtered_distance


def count_sectors(angle):
    # Count total recruitment values in sectors of 10 degrees bandwidth
    sectors = np.zeros(len(angle))
    for a in range(0, len(angle)):
        sectors[a] = int(angle[a] / 10)
    return sectors


def connectivity_by_sectors(sectors, connectivity, distance):
    # Calculate the total connectivity by sector
    sectors = np.array(sectors)
    connectivity = np.array(connectivity)
    distance = np.array(distance)
    vector = np.vectorize(np.int_)
    unique_sector = np.unique(sectors)
    sum_connectivity = []
    total_reefs = []
    eff_reefs = []
    avg_distance = []
    bandwidth = []
    weight_avg = []
    for band in unique_sector:
        index = np.array(np.argwhere(sectors == band))
        index = vector(index)
        total_reefs.append(len(index))
        sum_connectivity.append(connectivity[index].sum())
        if (np.sum(connectivity[index]) == 0):
            weights = np.zeros(len(index))
            eff_reefs.append(0)
        else:
            weights = (connectivity[index] /
                       np.sum(connectivity[index])).ravel()
            eff_reefs.append(np.count_nonzero(connectivity[index]))
        weight_avg.append(np.sum(distance[index].ravel() * weights))
        avg_distance.append(distance[index].mean())
        bandwidth.append(band)
    return np.array(sum_connectivity), np.array(bandwidth), np.array(avg_distance), np.array(weight_avg), np.array(eff_reefs), np.array(total_reefs)


def select_top_sectors(connectivity_array, num_sectors):
    # Select the num_sectors sectors with the highest connectivity value
    connectivity_array = np.array(connectivity_array)
    top_sectors = sorted(connectivity_array, reverse=True)[:num_sectors]
    top_sectors = np.trim_zeros(top_sectors)
    selected = np.where(np.isin(connectivity_array, top_sectors))[0]
    return np.array(selected)


def angle_circ_distance(angle_one, angle_two):
    """
    This function calculates the distance between to angles.
    """
    return (1-(np.cos(angle_one * math.pi / 180 - angle_two * math.pi / 180)))


def sort_angles(selected_bandwidth):
    """"
    This function arrange the vector of angles by distance between each other
    """
    selected_bandwidth = selected_bandwidth * 10
    sorted_angle_array = np.zeros(len(selected_bandwidth))
    sorted_angle_array[0] = selected_bandwidth[0]
    selected_bandwidth = np.delete(selected_bandwidth, 0)
    current = 0
    effective_distance = []
    while (len(selected_bandwidth) > 0):
        distances = []
        for i in range(0, len(selected_bandwidth)):
            distances.append(angle_circ_distance(
                sorted_angle_array[current], selected_bandwidth[i]))
        next_in_line = np.argmin(distances)
        effective_distance.append(np.min(distances))
        current += 1
        sorted_angle_array[current] = selected_bandwidth[next_in_line]
        selected_bandwidth = np.delete(selected_bandwidth, next_in_line)
    sorted_angle_array = [int(angle/10) for angle in sorted_angle_array]
    return (np.array(sorted_angle_array), np.array(effective_distance))


def select_boundary_sectors(top_sectors, total_connectivity):
    """"
    Select the boundary sectors around the top sectors
    """
    selected = np.array(
        np.unique(np.array((top_sectors + 1, top_sectors - 1)).ravel()))
    # Remove duplicated and negative index
    if (any(selected < 0)):
        selected[np.where(selected < 0)] = len(total_connectivity)-1
    if (any(selected >= len(total_connectivity))):
        selected[np.where(selected >= len(total_connectivity))] = 0
    selected = np.unique(selected)
    selected = selected[~np.in1d(selected, top_sectors)]
    temp_connectivity = sorted(
        total_connectivity[np.array(selected)], reverse=True)[0]
    final_select = np.array(
        np.where((temp_connectivity == total_connectivity)))
    return (final_select)


def find_two_groups(sectors_array):
    """
    Finds two groups of nearest values in a given list.
    """
    angles, distances = sort_angles(sectors_array)
    angle_sector_diff = np.diff(angles)
    # if the distance is higher than 1 (which means 10 degrees), it should find 2 groups.
    if (len(angles) > 1 and max(angle_sector_diff) > 1):
        cutting = int(np.max(np.where(distances == distances.max()))) + 1
        return angles[:cutting], angles[cutting:]
    else:
        return angles, []


def find_positions(list_of_sectors, sectors):
    """
    Finds the positions of specific values within a list.
    """
    if (len(sectors) > 1):
        first_sector = np.where(np.in1d(list_of_sectors, sectors[0]))[0]
        second_sector = np.where(np.in1d(list_of_sectors, sectors[1]))[0]
        return first_sector, second_sector
    else:
        only_sector = np.where(np.in1d(list_of_sectors, sectors[0]))[0]
        return only_sector, []


def calculate_ds(vector_angles):
    """
    Find the total angle from the source (DS) and the distance from the angle 0 (S)
    """
    DS = np.max(vector_angles) - np.min(vector_angles)
    if (DS == 0):
        DS = 10
    S = np.min(vector_angles) + (DS / 2)
    if (np.min(vector_angles) < np.max(vector_angles) - 180):
        DS = 360 - np.max(vector_angles) + np.min(vector_angles)
        S = (np.max(vector_angles) + (DS / 2)) % 360
    return (DS, S)


# Load input files for connectivity, angle, and shape information
# Step 1: Load shapefile and get centroids and reef order
shapefile = '../datasets/shapefiles/AIMS_shapefile/Simplified_LTMS_AIMS.shp'
centroids, reef_order = get_centroids(shapefile)

# Step 2: Load distance and angle matrices
distance_matrix = np.loadtxt('../outputs/GBR_zones_distance.csv',
                             delimiter=',', skiprows=1, usecols=range(1, len(reef_order)+1))
angle_matrix = np.loadtxt('../outputs/GBR_zones_angles.csv',
                          delimiter=',', skiprows=1, usecols=range(1, len(reef_order)+1))

distance_matrix_ordered = distance_matrix
angle_matrix_ordered = angle_matrix

# Step 3: Main loop over the connectivity files

nc_dataset = xr.open_dataset('../outputs/connectivity_matrices.nc')
## then the length of the time dimension is the number of time steps
for time_step in nc_dataset['time']:
    period = time_step.values.astype('datetime64[D]')
    ## just the year month and day
    connectivity_matrix = nc_dataset['connectivity'].sel(time=time_step).values
    #distance_matrix_ordered, angle_matrix_ordered = reef_ordering(
    #    distance_matrix, angle_matrix, reef_order)
    # Main code:
    # Step 4: Open output file for writing kernel parameters
    Kernel_outFile = open(str(period) + 'Kernel_parameters_Corals_sink.csv', 'w')
    Kernel_outFile.write("reef_ID,lat,lon,connectivity_sector_01,S_sector_01,DS_sector_01,Distance_sector_01,Proportion_reefs_01,connectivity_sector_02,S_sector_02,DS_sector_02,Distance_sector_02,Proportion_reefs_02\n")
    # Step 5: Iterate over sink areas
    for source_area in range(0, len(connectivity_matrix[0])):
        Kernel_outFile.write(str(
            source_area) + ',' + str(centroids[source_area].y) + ',' + str(centroids[source_area].x) + ',')
        # Step 6: Select source locations based on connectivity, angle, and distance matrices
        source_connect, source_angle, source_distance = select_sink_locations(
            connectivity_matrix, angle_matrix_ordered, distance_matrix_ordered, source_area)
        # source_connect, source_angle, source_distance = select_source_locations(
        #     connectivity_matrix, angle_matrix_ordered, distance_matrix_ordered, source_area)
        # Step 7: Remove outliers from the selected sources
        filter_connect, filter_angle, filter_centroid, filter_distance = remove_outliers(
            source_connect, source_angle, centroids, source_distance, percent=5)
        # Step 8: Count sectors based on filtered angle values
        sectors_array = count_sectors(filter_angle)
        # Step 9: Calculate total connectivity, bandwidth, average distance, and weighted distance by sectors
        total_connectivity, bandwidth, avg_distance, wgt_distance, effective_reefs, total_reefs = connectivity_by_sectors(
            sectors_array, filter_connect, filter_distance)
        # Step 10: Check if total connectivity is zero, write 'nan' values, and continue to the next sink area
        if (sum(total_connectivity) == 0):
            Kernel_outFile.write('nan,nan,nan,nan,nan,nan,nan,nan,nan,nan\n')
            continue
        n_sectors = 2  # number of sectors
        # Step 11: Select top sectors based on total connectivity
        selected_sector = select_top_sectors(total_connectivity, n_sectors)
        # Step 12: Calculate proportion of recruitment and add additional sectors if needed
        prop_recruitment = sum(
            total_connectivity[selected_sector])/sum(total_connectivity)
        while prop_recruitment < 0.8:
            selected = select_boundary_sectors(
                selected_sector, total_connectivity)
            selected_sector = np.append(selected_sector, selected)
            prop_recruitment = sum(
                total_connectivity[selected_sector])/sum(total_connectivity)
        # 12.5 cleaning the vector with zeros
        selected_sector = selected_sector[np.where(
            total_connectivity[selected_sector])[0]]
        # Step 13: Find two groups based on bandwidth of selected sectors
        bandwidth_sectors = find_two_groups(bandwidth[selected_sector])
        bandwidth_sectors_filtered = list(
            filter(lambda x: len(x) > 0, bandwidth_sectors))
        # Step 14: Find positions in the bandwidth matrix for sector 1 and sector 2
        sector_1, sector_2 = find_positions(
            bandwidth, bandwidth_sectors_filtered)
        # Step 15: Write parameters for sector 1
        ds1, s1 = calculate_ds(bandwidth[sector_1] * 10)
        Kernel_outFile.write(str(sum(total_connectivity[sector_1])) + ',' + str(s1) + ',' + str(ds1) + ',' + str(
            np.mean(wgt_distance[sector_1])) + ',' + str(np.sum(effective_reefs[sector_1]) / np.sum(total_reefs[sector_1])) + ',')
        # Step 16: Write parameters for sector 2 if it exists
        if (len(sector_2) > 0):
            ds2, s2 = calculate_ds(bandwidth[sector_2] * 10)
            Kernel_outFile.write(str(sum(total_connectivity[sector_2])) + ',' + str(s2) + ',' + str(ds2) + ',' + str(
                np.mean(wgt_distance[sector_2])) + ',' + str(np.sum(effective_reefs[sector_2]) / np.sum(total_reefs[sector_2])) + '\n')
        else:
            Kernel_outFile.write('nan,nan,nan,nan,nan\n')
    Kernel_outFile.close()
