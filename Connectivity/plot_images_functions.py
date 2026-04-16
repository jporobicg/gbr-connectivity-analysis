import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_probability_of_connection(nc_dataset, gdf, output_file):
    # Sum connectivity values across a dimension, e.g., 'target'
    connectivity_sum = nc_dataset['connectivity'].sum(dim='sink')

    # Convert sums to probabilities
    total_connectivity = connectivity_sum.sum()
    probabilities = connectivity_sum / total_connectivity

    # Convert to pandas DataFrame for plotting
    prob_df = probabilities.to_dataframe(name='probability').reset_index()
    # calculate the mean probability of connection by source
    mean_probabilities = prob_df.groupby('source')['probability'].mean().reset_index()
    # using this probablilities we can plot the probability of connection by source using the shapefile
    # Ensure 'FID' is of the same type as 'source' in mean_probabilities for a successful merge
    gdf = gdf.copy()
    gdf['FID'] = gdf['FID'].astype(str)
    mean_probabilities['source'] = mean_probabilities['source'].astype(str)

    # Merge the GeoDataFrame with mean_probabilities DataFrame
    merged_gdf = gdf.merge(mean_probabilities, left_on='FID', right_on='source')

    # Generate a color palette: normalize probabilities to [0, 1]
    prob_normalized = (mean_probabilities['probability'] - mean_probabilities['probability'].min()) / (mean_probabilities['probability'].max() - mean_probabilities['probability'].min())
    cmap = sns.color_palette("inferno_r", as_cmap=True)
    colors = list(cmap(prob_normalized))

    fig, axs = plt.subplots(1, 2, figsize=(12, 8))

    # Bar plot of mean probabilities
    sns.barplot(x='source', y='probability', data=mean_probabilities, hue='source', palette=dict(zip(mean_probabilities['source'], colors)), legend=False, ax=axs[0])
    axs[0].set_title('Probability of Connection', fontsize=16)
    axs[0].set_xlabel('Source', fontsize=14)
    axs[0].set_ylabel('Mean Probability', fontsize=14)
    plt.xticks(rotation=45)

    # Spatial distribution
    merged_gdf.plot(column='probability', cmap='inferno_r', legend=True, ax=axs[1])
    axs[1].set_title('Probability of Connection', fontsize=16)
    merged_gdf.boundary.plot(ax=axs[1], edgecolor='black', linewidth=1)
    axs[1].set_xlabel('Longitude', fontsize=14)
    axs[1].set_ylabel('Latitude', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def heatmap_connectivity_by_year(nc_dataset, output_file):
    # Time has 3 events per year: sum connectivity within each year first
    conn = nc_dataset['connectivity']
    years = conn.time.dt.year
    connectivity_by_year = conn.groupby(years).sum(dim='time')  # (year, source, sink)
    years_list = connectivity_by_year.coords['year'].values
    n_years = len(years_list)
    n_cols = min(3, n_years)
    n_rows = (n_years + n_cols - 1) // n_cols
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_years == 1:
        axs = np.array([axs])
    axs = axs.flatten()
    for i, year in enumerate(years_list):
        mat = connectivity_by_year.sel(year=int(year)).values.squeeze()  # 2D (source, sink)
        sns.heatmap(mat, ax=axs[i], cmap='inferno_r', cbar_kws={'label': 'Connectivity'})
        axs[i].set_title(f'Year {int(year)}')
    for j in range(i + 1, len(axs)):
        axs[j].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def heatmap_connectivity_by_year_by_region(
    nc_dataset,
    gdf,
    output_file,
    mapping_path=None,
    coral_shapefile_path=None,
):
    """
    Heatmap of connectivity summed by year, aggregated by region (FID).
    Reef-level connectivity (3806x3806) is aggregated to region-level (30x30)
    using the mapping: coral reef index -> ltms_id (region FID).
    """
    import geopandas as gpd
    from pathlib import Path

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    if mapping_path is None:
        mapping_path = project_root / "datasets/reefs/combined_mapping_index.txt"
    if coral_shapefile_path is None:
        coral_shapefile_path = project_root / "datasets/shapefiles/gbr1_coral_1m_merged_buffer0p001.shp"

    coral_gdf = gpd.read_file(coral_shapefile_path)
    mapping = pd.read_csv(mapping_path, sep="\t")
    coral_id_to_ltms = mapping.set_index("coral_id")["ltms_id"].to_dict()
    index_to_region = np.array(
        [coral_id_to_ltms.get(coral_gdf["FID"].iloc[i], -1) for i in range(coral_gdf.shape[0])]
    )
    region_fids = sorted(gdf["FID"].unique())
    n_regions = len(region_fids)
    fid_to_idx = {int(fid): k for k, fid in enumerate(region_fids)}

    # Map each reef index to region index (-1 if unmapped)
    src_region_idx = np.array([fid_to_idx.get(int(r), -1) for r in index_to_region])
    sink_region_idx = src_region_idx.copy()

    conn = nc_dataset["connectivity"]
    years = conn.time.dt.year
    connectivity_by_year = conn.groupby(years).sum(dim="time")
    years_list = connectivity_by_year.coords["year"].values
    n_years = len(years_list)

    conn_by_region = np.zeros((n_years, n_regions, n_regions))
    for t in range(n_years):
        mat = connectivity_by_year.sel(year=int(years_list[t])).values.squeeze()
        valid = (src_region_idx >= 0) & (sink_region_idx >= 0)
        src_2d = np.broadcast_to(src_region_idx[:, None], mat.shape)
        sink_2d = np.broadcast_to(sink_region_idx[None, :], mat.shape)
        valid_2d = (src_2d >= 0) & (sink_2d >= 0)
        np.add.at(
            conn_by_region[t],
            (src_2d[valid_2d], sink_2d[valid_2d]),
            mat[valid_2d],
        )

    n_cols = min(3, n_years)
    n_rows = (n_years + n_cols - 1) // n_cols
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_years == 1:
        axs = np.array([axs])
    axs = axs.flatten()
    for i in range(n_years):
        sns.heatmap(
            conn_by_region[i],
            ax=axs[i],
            cmap="inferno_r",
            cbar_kws={"label": "Connectivity"},
            xticklabels=region_fids,
            yticklabels=region_fids,
        )
        axs[i].set_title(f"Year {int(years_list[i])}")
        axs[i].set_xlabel("Sink region (FID)")
        axs[i].set_ylabel("Source region (FID)")
    for j in range(n_years, len(axs)):
        axs[j].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

