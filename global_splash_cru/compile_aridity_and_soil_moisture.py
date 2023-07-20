from pathlib import Path
import argparse

import xarray
import pandas
import numpy as np
from pyrealm.pmodel.functions import calc_soilmstress_stocker, calc_soilmstress_mengoli


def compile_aridity_and_soil_moisture(basename: Path):
    """Import precipitation and SPLASH PET for a set of sites and reduce to annual
    estimates of the aridity index."""

    # Get a list of the output files and open them
    results_files = Path(basename / "splash_outputs").glob("results_block_*.nc")
    data = xarray.open_mfdataset(results_files)

    # Sort the dataset  to match the grid order
    data = data.sortby("cell_id")

    # Load WFDE5 elevation to add lat and lon data indexed by cell ids
    elev = xarray.load_dataarray(basename.parent / "data" / "Asurf_land_cells.nc")
    elev = elev.where(elev.cell_id.isin(data.cell_id), drop=True)

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

    # ----------------------------------
    # Export gridded aridity index
    # ----------------------------------

    # Annual total precipitation and PET and hence aridity index by year and across
    # climatology - setting np.inf to np.nan
    annual_AI_sum_meth = (
        data["pet"].groupby("time.year").sum() / data["pre"].groupby("time.year").sum()
    )
    annual_AI_sum_meth = annual_AI_sum_meth.where(annual_AI_sum_meth < np.inf, np.nan)

    # 20 year climatological
    clim_AI_sum_meth = data["pet"].sum(dim="time") / data["pre"].sum(dim="time")
    clim_AI_sum_meth = clim_AI_sum_meth.where(clim_AI_sum_meth < np.inf, np.nan)

    aridity = xarray.Dataset(
        {"annual_ai": annual_AI_sum_meth, "climatological_ai": clim_AI_sum_meth}
    )

    # Export AI summary
    aridity.to_netcdf(basename / "aridity_index_data.nc")

    # ----------------------------------
    # Soil moisture and soilmstress to grid
    # ----------------------------------

    # Calculate soil moisture stress - broadcast ai data to same shape as wn
    # and calculate the relative soil moisture using the SPLASH fixed bucket size.
    ai = clim_AI_sum_meth.broadcast_like(data["wn"]).compute().data
    soilm = data["wn"].compute()
    relative_soilm = soilm.data / 150
    soilmstress_mengoli = calc_soilmstress_mengoli(
        aridity_index=ai, soilm=relative_soilm
    )

    soilm_data = xarray.Dataset(
        {
            "soilm": soilm,
            "soilmstress_mengoli": xarray.DataArray(
                soilmstress_mengoli, dims=soilm.dims, coords=soilm.coords
            ),
        }
    )

    # Export soil moisture data by year
    years, datasets = zip(*soilm_data.groupby("time.year"))
    paths = [basename / "soil_moisture_grids" / f"soil_moisture_{y}.nc" for y in years]
    xarray.save_mfdataset(datasets, paths)


if __name__ == "__main__":
    # Parse the arguments provided to main and execute the function.
    parser = argparse.ArgumentParser("compile_aridity_and_soil_moisture")
    parser.add_argument("--basename", type=Path)
    args = parser.parse_args()

    compile_aridity_and_soil_moisture(basename=args.basename)
