import xarray
import numpy as np

vars = {
    "tmp": (
        "../data/cru_ts4.06.2001.2010.tmp.dat.nc",
        "../data/cru_ts4.06.2011.2020.tmp.dat.nc",
    ),
    "cld": (
        "../data/cru_ts4.06.01.2001.2010.cld.dat.nc",
        "../data/cru_ts4.06.01.2011.2020.cld.dat.nc",
    ),
    "pre": (
        "../data/cru_ts4.06.2001.2010.pre.dat.nc",
        "../data/cru_ts4.06.2011.2020.pre.dat.nc",
    ),
}


# Create a single integrated file for each variable.
for v_name, files in vars.items():
    # Open the datasets into a MF dataset and extract the main variable, dropping stn
    data = xarray.open_mfdataset(files)[v_name]

    # Chunk by time series to speed iteration over time series
    data.to_netcdf(f"../data/cru_2001_2020_{v_name}.nc")

# Load WFDE5 elevation and add a cell_id index across all cells
elev = xarray.load_dataarray("../data/ASurf_WFDE5_CRU_v2.0.nc")
elev = elev.stack(cell_id=("lat", "lon"))

# Convert multi-index back an assign individual cell ids
elev = elev.reset_index("cell_id")
elev = elev.assign_coords({"cell_id": np.arange(np.product(elev.shape))})

# Reduce to only land cells and save
elev = elev.loc[np.logical_not(np.isnan(elev.values))]
elev.to_netcdf("../data/Asurf_land_cells.nc")
