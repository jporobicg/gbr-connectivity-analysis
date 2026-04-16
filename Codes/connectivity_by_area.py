## calculate new connectivity matrices By Management areas (30x30)
import xarray as xr
import pandas as pd
import numpy as np
import glob
import re
from datetime import datetime

## Functions
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
    new_connectivity_matrix = pd.DataFrame(0.0, index=unique_values, columns=unique_values, dtype=float)
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


def create_mapped_connectivity_matrix_from_netcdf(
    nc_path,
    mapping_index_path,
    output_path=None,
    connectivity_var="connectivity",
    time_dim="time",
    slicing=False,
):
    """
    Read a NetCDF with connectivity matrices (time, source, target/sink), apply the
    same mapping as create_mapped_connectivity_matrix to each time step, and return
    (and optionally save) a new xarray Dataset with mapped dimensions (ltms_id).

    Parameters:
    - nc_path (str or path-like): Path to the input NetCDF file.
    - mapping_index_path (str or path-like): Path to the mapping index CSV/TSV
      (columns: coral_id, ltms_id).
    - output_path (str or None): If set, save the mapped dataset to this path.
    - connectivity_var (str): Name of the connectivity variable in the NetCDF.
    - time_dim (str): Name of the time dimension.

    Returns:
    - xr.Dataset: Mapped connectivity dataset with dims (time, source, target).
    """
    mapping_index = pd.read_csv(mapping_index_path, sep="\t")
    ds = xr.open_dataset(nc_path)
    if slicing is not False:
        ds = ds.isel(treatment=0,sample=0)

    conn = ds[connectivity_var]

    # Support both 'target' and 'sink' dimension names
    sink_dim = "sink"
    source_dim = "source"

    times = conn.coords[time_dim].values
    source_ids = conn.coords[source_dim].values
    sink_ids = conn.coords[sink_dim].values

    data_arrays = []
    for t in times:
        mat = conn.sel({time_dim: t}).values  # 2D: (source, sink)
        connectivity_matrix = pd.DataFrame(
            mat, index=source_ids, columns=sink_ids
        )
        new_connectivity_matrix = create_mapped_connectivity_matrix(
            connectivity_matrix, mapping_index
        )
        data_arrays.append(new_connectivity_matrix.values)

    new_source = np.array(new_connectivity_matrix.index)
    new_sink = np.array(new_connectivity_matrix.columns)

    ds_mapped = xr.Dataset(
        {
            "connectivity": (
                (time_dim, "source", "sink"),
                np.array(data_arrays),
            )
        },
        coords={
            time_dim: times,
            "source": new_source,
            "sink": new_sink,
        },
    )

    # Preserve/copy metadata
    src_conn = ds[connectivity_var]
    ds_mapped["connectivity"].attrs["description"] = (
        src_conn.attrs.get("description", "Connectivity matrix by date (mapped to LTMS areas)")
    )
    ds_mapped["connectivity"].attrs["units"] = (
        src_conn.attrs.get("units", "number of individuals")
    )
    for key in ["title", "institution", "source", "history", "shapefiles"]:
        if key in ds.attrs:
            ds_mapped.attrs[key] = ds.attrs[key]
    ds_mapped.attrs["history"] = ds.attrs.get("history", "") + "; Mapped by create_mapped_connectivity_matrix_from_netcdf"

    ds.close()
    if output_path:
        ds_mapped.to_netcdf(output_path)
    return ds_mapped


nc_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/connectivity_matrices/connectivity_merulinidae_single.nc'
mapping_index_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/reefs/combined_mapping_index.txt'
output_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/connectivity_matrices/connectivity_merulinidae_single_by_area.nc'
create_mapped_connectivity_matrix_from_netcdf(nc_path, mapping_index_path, output_path, slicing=True)