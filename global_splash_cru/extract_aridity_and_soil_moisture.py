from pathlib import Path
import argparse

import xarray
import pandas
import numpy as np


def extract_annual_aridity(basename: Path):
    """Import precipitation and SPLASH PET for a set of sites and reduce to annual
    estimates of the aridity index."""

    # Get a list of the output files and a dict to store summaries
    results_files = basename / Path("results").glob("results_block_*.nc")
    sections = {}
    soilm_sections = {}

    for each_file in results_files:
        # Open the block file
        print(each_file)
        data = xarray.open_dataset(each_file)

        # Get the data for precipitation, PET and soil moisture
        precip = data["pre"]
        pet = data["pet"]

        # Annual total precipitation and PET and hence aridity index by year
        annual_AI_sum_meth = (
            pet.groupby("time.year").sum() / precip.groupby("time.year").sum()
        )

        # 20 year climatological
        clim_AI_sum_meth = pet.sum(dim="time") / precip.sum(dim="time")

        # # Calculations using daily average values
        # daily_AI = pet / precip
        # annual_AI_mean_meth = daily_AI.groupby("time.year").mean()
        # clim_AI_mean_meth = daily_AI.mean(dim="time")

        sections[each_file] = xarray.Dataset(
            {
                "annual_AI_sum_meth": annual_AI_sum_meth,
                # "annual_AI_mean_meth": annual_AI_mean_meth,
                "clim_AI_sum_meth": clim_AI_sum_meth,
                # "clim_AI_mean_meth": clim_AI_mean_meth,
            }
        )

        soilm_sections[each_file] = data["wn"]
        data.close()

    # ----------------------------------
    # Aridity index to grid
    # ----------------------------------

    # Compile the sections into a single dataset along the cell_id axis and then
    # calculate climatological AI across years
    data = xarray.concat(sections.values(), dim="cell_id")

    # Sort the dataset  to match the grid order
    data = data.sortby("cell_id")

    # Load WFDE5 elevation to add lat and lon data indexed by cell ids
    elev = xarray.load_dataarray("Asurf_land_cells.nc")
    data["lat"] = xarray.DataArray(elev.lat.data, dims=("cell_id",))
    data["lon"] = xarray.DataArray(elev.lon.data, dims=("cell_id",))

    # Unstack lat and lon back to 2 dimensions from cell id
    # https://stackoverflow.com/questions/70861487/
    data = (
        data.assign_coords(
            {
                "cell_id": pandas.MultiIndex.from_arrays(
                    [data.lat.data, data.lon.data],
                    names=["lat", "lon"],
                )
            }
        )
        .drop_vars(["lat", "lon"])
        .unstack("cell_id")
    )

    # insert missing latitudes
    # https://stackoverflow.com/questions/68207994
    all_lats = {"lat": np.arange(-89.75, 90, 0.5)}
    data = data.reindex(all_lats, fill_value=np.nan)

    # Set odd infinity values to np.nan
    clim_data = data["clim_AI_sum_meth"]
    data["clim_AI_sum_meth"] = clim_data.where(clim_data < np.inf, np.nan)
    data["log_clim_AI_sum_meth"] = data["clim_AI_sum_meth"].log()

    # Export AI summary
    data.to_netcdf(basename / "aridity_index_data.nc")

    # ----------------------------------
    # Soil moisture to grid
    # ----------------------------------

    # Compile the sections into a single dataset along the cell_id axis and then
    # calculate climatological AI across years
    soilm_data = xarray.concat(soilm_sections.values(), dim="cell_id")

    # Sort the dataset  to match the grid order
    soilm_data = soilm_data.sortby("cell_id")

    # Load WFDE5 elevation to add lat and lon data indexed by cell ids
    soilm_data["lat"] = xarray.DataArray(elev.lat.data, dims=("cell_id",))
    soilm_data["lon"] = xarray.DataArray(elev.lon.data, dims=("cell_id",))

    # Unstack lat and lon back to 2 dimensions from cell id
    # https://stackoverflow.com/questions/70861487/
    soilm_data = (
        soilm_data.assign_coords(
            {
                "cell_id": pandas.MultiIndex.from_arrays(
                    [soilm_data.lat.data, soilm_data.lon.data],
                    names=["lat", "lon"],
                )
            }
        )
        .drop_vars(["lat", "lon"])
        .unstack("cell_id")
    )

    # insert missing latitudes
    # https://stackoverflow.com/questions/68207994
    all_lats = {"lat": np.arange(-89.75, 90, 0.5)}
    soilm_data = soilm_data.reindex(all_lats, fill_value=np.nan)

    # Export soil moisture data
    soilm_data.to_netcdf(basename / "soilm_data.nc")


if __name__ == "__main__":
    # Parse the arguments provided to main and execute the function.
    parser = argparse.ArgumentParser("extract_annual_aridity")
    parser.add_argument("--basename", type=Path)
    args = parser.parse_args()

    extract_annual_aridity(basename=args.basename)
