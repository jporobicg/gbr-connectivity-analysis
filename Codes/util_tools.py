import xarray as xr
import pandas as pd
import numpy as np
from scipy.integrate import quad
import sys
from shapely.geometry import Point
import shapely
import math
from tqdm import tqdm
import pandas as pd
from numba import jit, njit
# Numba translates Python functions to optimized machine code at runtime.
import numba
import os
from joblib import Parallel, delayed, parallel_backend
import geopandas as gpd
from shapely.geometry import LineString

## ~~~~~~~~~~~~~~~~~~ ##
## ~     Functions  ~ ##
## ~~~~~~~~~~~~~~~~~~ ##


def bathtub_curve(lmbda, v, sigma):
    """
    Transforms the bathtub equation into a function that can be integrated
    using the trapezoidal function.

    Parameters
    ----------
    lmbda : float
        The scale parameter of the Weibull distribution used in the bathtub equation.
        Must be a positive integer or float.
    v : float
        The shape parameter of the Weibull distribution used in the bathtub equation.
        Must be a positive integer or float.
    sigma : float
        A constant that determines the shape of the bathtub curve. Must be a positive
        integer or float. It can be zero, becoming an exponential function.

    Returns
    -------
    function
        A lambda function that represents the bathtub curve equation. The function takes
        a single argument, t, which represents the age of the coral. The function returns
        the value of the bathtub curve equation at that age.
    """

    def u(t): return (lmbda * v * pow((lmbda * t), v - 1)) / \
        (1 - sigma * pow((lmbda * t), v))
    return (u)


def piecewise_decay(ages, Tcp, lmbda1, lmbda2, v1, v2, sigma1, sigma2):
    """
    Calculates the probability of survival of corals larvaes at different ages, using the piecewise
    Weibull-Weibull survival model described by Moneghetti et al. (2019).

    Parameters
    ----------
    ages : list
        A list of ages of the corals. Each age must be a positive integer or float.
    Tcp : float
        The age (inflection point) at which the corals transition from a Weibull survival curve to another Weibull
        survival curve. Must be greater than 0 and less than the maximum age in `ages`.
    lmbda1 : float
        The scale parameter of the Weibull survival curve in the first phase. Must be a positive
        integer or float.
    lmbda2 : float
        The scale parameter of the Weibull survival curve in the second phase. Must be a positive
        integer or float.
    v1 : float
        The shape parameter of the Weibull survival curve in the first phase. Must be a positive
        integer or float.
    v2 : float
        The shape parameter of the Weibull survival curve in the second phase. Must be a positive
        integer or float.
    sigma1 : float
        The standard deviation of the Gaussian noise added to the survival curve in the first phase.
        Must be a positive integer or float. It can be zero, becoming an exponential function.
    sigma2 : float
        The standard deviation of the Gaussian noise added to the survival curve in the second phase.
        Must be a positive integer or float. It can be zero, becoming an exponential function.

    Returns
    -------
    list
        A list of survival probabilities for the corals, calculated using the piecewise
        Weibull-Weibull survival model. Each survival probability corresponds to the age
        in the input list `ages`.
    """
    fx1 = bathtub_curve(lmbda1, v1, sigma1)
    fx2 = bathtub_curve(lmbda2, v2, sigma2)
    decay = []
    for age in range(0, len(ages)):
        if (ages[age] < Tcp):
            area = quad(fx1, 0, ages[age])[0]
        else:
            area = quad(fx1, 0, Tcp)[0] + quad(fx2, Tcp, ages[age])[0]
        decay.append(math.exp(-area))
    return decay


def exponential_growth_function(initial_value, growth_rate, time_period):
    """
    Calculates the exponential growth over a given time period.

    Parameters:
    - initial_value (float): The starting value of the exponential growth.
    - growth_rate (float): The rate at which the value is growing, expressed as a decimal or fraction.
    - time_period (int): The number of time periods over which the growth occurs (larval age in days).

    Returns:
    - float: The value after the exponential growth.
    """
    return initial_value * pow((1 + growth_rate), time_period)


def exponential_decay_function(initial_value, decay_rate, time, initial_time=0.0):
    """
    Calculate the value of an exponentially decaying quantity at a given time.
    This is a very simple decay function assuming mean lifetime decay for the entire particle
    Parameters:
    - initial_value: float
        The initial value of the decaying quantity (this is set to 1 as is a proportion).
    - decay_rate: float
        The rate at which the quantity decays per unit time.
    - time: float
        The time at which the decay is evaluated (age of the larvae).
    Returns:
    - result: float
       The value of the decaying quantity at the given time.
    """
    return initial_value * np.exp(-decay_rate * (time - initial_time))


def exponential_competence(ages, tc, peak_age):
    """
    Calculates the competence level based on age using an exponential growth function.

    Parameters:
    - ages (list): A list of ages.
    - tc (float): The threshold age below which the competence level is 0.
    - peak_age (float): The age at which the competence level reaches its peak (1).

    Returns:
    - list: A list of competence levels corresponding to each age.
    """
    competence = []
    for age in range(0, len(ages)):
        if (ages[age] < tc):
            probability = 0
        elif (ages[age] < peak_age):
            probability = exponential_growth_function(0.01, 0.245, ages[age])
        else:
            probability = 1
        competence.append(probability)
    return competence


def exponential_decay(ages, tc, phases_age):
    """
    Calculates the decay probability based on age using an exponential decay function.

    Parameters:
    - ages (list): A list of ages.
    - tc (float): The threshold age below which the decay probability is 0.
    - phases_age (float): The age at which the decay probability transitions to a different rate.

    Returns:
    - list: A list of decay probabilities corresponding to each age.

    """
    decay = []
    for age in range(0, len(ages)):
        if (ages[age] < tc):
            probability = 0
        elif (ages[age] < phases_age):
            probability = exponential_decay_function(1.0, 0.2, ages[age])
        else:
            probability = exponential_decay_function(1.0, 0.5, ages[age])
        decay.append(probability)
    return decay


def piecewise_competence(ages, tc, Tcp, alpha, beta1, beta2, v):
    """
    Calculates the larval competence values at different ages (days), using the piecewise
    Weibull-exponential competence model. This function is a replica of the R code used by
    Moneghetti et al. (2019) to calculate competence.

    Parameters
    ----------
    ages : list
        A list of larvaes ages. Each age must be a positive integer or float.
    tc : float
        The age at which the larvaes reaches their maximum competence level. Must be a
        positive integer or float.
    Tcp : float
        The age at which the larvaes starts to experience a decline in competence. Must be
        greater than tc and a positive integer or float.
    alpha : float
        The scale parameter of the Weibull distribution. Must be a positive integer or float.
    beta1 : float
        The shape parameter of the Weibull distribution in the early decline phase. Must be
        a positive integer or float.
    beta2 : float
        The shape parameter of the Weibull distribution in the late decline phase. Must be
        a positive integer or float.
    v : float
        The exponential decay parameter in the early decline phase. Must be a positive
        integer or float.

    Returns
    -------
    list
        A list of competence values for larvaes, calculated using the piecewise
        Weibull-exponential competence model. Each competence value corresponds to the age
        in the input list `ages`.
    """
    competence = []
    for age in range(0, len(ages)):
        if (ages[age] < tc):
            area = 0
        if (ages[age] >= tc and ages[age] <= Tcp):
            def fxtau_early(tau): return alpha * math.exp(-alpha *
                                                          (tau - tc)) * math.exp(- pow((beta1 * (ages[age]-tau)), v))
            area = quad(fxtau_early, tc, ages[age])[0]
        if (ages[age] > Tcp):
            def fxtau_late_first(tau): return alpha * math.exp(-alpha * (tau - tc)) * \
                math.exp(- pow((beta1 * (Tcp-tau)), v)) * \
                math.exp(-beta2 * (ages[age] - Tcp))
            def fxtau_late_second(tau): return alpha * math.exp(-alpha *
                                                                (tau - tc)) * math.exp(- beta2 * (ages[age]-tau))
            area = quad(fxtau_late_first, tc, Tcp)[
                0] + quad(fxtau_late_second, Tcp, ages[age])[0]
        competence.append(area)
    return (competence)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ In polygon algorithm and optimizers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


@jit(nopython=True)
def point_in_polygon(x, y, polygon):
    """
    Determines whether a point is inside a polygon using the ray-casting algorithm.

    Parameters
    ----------
    x : float
        The x-coordinate of the point to test.
    y : float
        The y-coordinate of the point to test.
    polygon : list of tuples
        A list of tuples representing the vertices of the polygon, in order.

    Returns
    -------
    bool
        True if the point is inside the polygon, False otherwise.
    """
    num_vertices = len(polygon)
    is_inside = False
    previous_x, previous_y = polygon[0]
    # start outside the polygon
    intersection_x = 0.0
    for i in numba.prange(num_vertices + 1):
        current_x, current_y = polygon[i % num_vertices]
        if y > min(previous_y, current_y):
            if y <= max(previous_y, current_y):
                if x <= max(previous_x, current_x):
                    if previous_y != current_y:
                        intersection_x = (
                            y - previous_y) * (current_x - previous_x) / (current_y - previous_y) + previous_x
                    if previous_x == current_x or x <= intersection_x:
                        is_inside = not is_inside
        previous_x, previous_y = current_x, current_y

    return is_inside


@njit(parallel=False)
def points_in_polygon(xs, ys, miny, maxy, polygon):
    """
    Test whether a point is inside a given polygon using the point-in-polygon algorithm.

    This function tests each point in the `points` array to determine if it is inside the polygon
    defined by the vertices in the `polygon` array. The function uses the point-in-polygon Ray casting
    algorithm to perform this test. The algorithm performs the even-odd-rule algorithm to find out
    whether a point is in a given polygon. This runs in O(n) time where n is the number of edges of the polygon.

    Parameters:
    -----------
    points: numpy.ndarray of shape (n, 2)
        Array of n points to test. Each point is defined by its x and y coordinates in columns 0 and 1 respectively.
    polygon: numpy.ndarray of shape (m, 2)
        Array of m vertices defining the polygon. Each vertex is defined by its x and y coordinates in columns 0 and 1 respectively.

    Returns:
    --------
    D: numpy.ndarray of shape (n,)
        Boolean array indicating whether each point in `points` is inside the polygon (`True`) or not (`False`).
    """
    D = np.empty(len(ys), dtype=numba.boolean)
    for i in range(len(D)):
        if ys[i] >= miny and ys[i] <= maxy:
            D[i] = point_in_polygon(xs[i], ys[i], polygon)
        else:
            D[i] = False
    return D

def find_missing_features(gdf_centroids, joined_gdf):
    """Find features in 'gdf_centroids' that are not in 'joined_gdf'."""
    missing_indices = set(gdf_centroids.index) - set(joined_gdf.index)
    print(f"Number of missing features: {len(missing_indices)}")
    # Convert missing_indices to a list before using it as an indexer
    missing_features = gdf_centroids.loc[list(missing_indices)]
    missing_features['closest_polygon_index'] = None
    return missing_features

def assign_closest_polygon(missing_features, GBR_gdf):
    """Assign the index of the closest polygon from 'GBR_gdf' to each missing feature."""
    for index, missing_feature in missing_features.iterrows():
        distances = GBR_gdf.geometry.distance(missing_feature.geometry)
        closest_polygon_index = distances.idxmin()
        missing_features.at[index, 'closest_polygon_index'] = GBR_gdf.iloc[closest_polygon_index]['FID']



def create_mapped_connectivity_matrix(connectivity_matrix, mapping_index):
    """
    Creates a new connectivity matrix based on a mapping dictionary, with unique values as indices and columns.
    The new matrix values are the average of the original matrix values mapped from the original to the new identifiers.

    Parameters:
    - connectivity_matrix (pd.DataFrame): Original connectivity matrix.
    - mapping_dict (dict): Mapping from original identifiers to new identifiers.

    Returns:
    - pd.DataFrame: New connectivity matrix with mapped identifiers and averaged values.
    """
    mapping_dict = mapping_index.set_index('coral_id')['ltms_id'].to_dict()
    # Extract unique values from mapping_dict.values() and sort them if needed
    unique_values = sorted(set(mapping_dict.values()))
    # Initialize new connectivity matrix with unique values as indices and columns
    new_connectivity_matrix = pd.DataFrame(index=unique_values, columns=unique_values).fillna(0)
    # Fill the new connectivity matrix based on the mapping_dict
    for row_new_id in unique_values:
        row_original_id = [key for key, value in mapping_dict.items() if value == row_new_id]
        for col_new_id in unique_values:
            col_original_id = [key for key, value in mapping_dict.items() if value == col_new_id]
            # Extract the original connectivity values
            values = connectivity_matrix.loc[row_original_id, col_original_id].values.flatten()
            # Calculate the average value, ignoring NaNs
            avg_value = np.nanmean(values) if values.size else 0
            # Fill the new connectivity matrix
            new_connectivity_matrix.loc[row_new_id, col_new_id] = avg_value
    return new_connectivity_matrix