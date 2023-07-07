"""Fitting the subdaily P Model to global data.

This script is to be used as part of an array job on the HPC to iterate over blocks of
cells from WFDE5 data in order to fit the subdaily P Model at global scale and 0.5°
resolution.
"""

import os
from pathlib import Path

import xarray
import numpy as np
import pandas as pd

from pyrealm.hygro import convert_sh_to_vpd
from pyrealm.pmodel.functions import calc_patm
from pyrealm.pmodel import PModelEnvironment, FastSlowPModel, FastSlowScaler


# ----------------------------------
# DATA PATHS
# ----------------------------------

# WFDE5 provides TEMP, PATM, VPD, PPFD
wfde_path = Path("/rds/general/project/lemontree/live/source/wfde5/wfde5_v2")
# SNU FAPAR downsampled from 0.05° to 0.5° gives FAPAR
fapar_path = Path("/rds/general/project/lemontree/live/derived/SNU_Ryu/FPAR_05d")
# NOAA global sea level time series provides CO2
co2_path = Path("/rds/general/project/lemontree/live/source/NOAA_CO2/co2_mm_gl.csv")


start_time = np.datetime64("2000-01-01")
end_time = np.datetime64("2020-01-01")

# ----------------------------------
# ELEVATION DATA
# ----------------------------------

# Read in the elevation data
asurf_data = xarray.load_dataarray(wfde_path / "Elev" / "ASurf_WFDE5_CRU_v2.0.nc")

# Use longitudinal stripes to control for timing: data is sampled at UTC, so local time
# varies with longitude and needs to be corrected. The blocksize sets the number of
# those vertical stripes to be loaded
stripe_width = 1
stripe_starts = np.arange(0, asurf_data.sizes["lon"], stripe_width)

# Select the stripe start for this subjob
array_index = os.environ["PBS_ARRAY_INDEX"]
lon_idx = stripe_starts[int(array_index)]

# https://stackoverflow.com/questions/66789660

# The WFDE5 v2 files are held in monthly files, so open all of the files as a dask
# multifile dataset and then access a spatial chunk of the data. It might be possible to
# improve the handling time by tuning the chunk sizes.
#
# The maximum time chunk size in any _one_ file is 31 days * 24 hours = 744 and the
# simulation currently needs to loop over longitudinal slices to apply appropriate
# offsets from UTC for local time.
#
# One longitudinal slice is (259400 * 360 * 1 * 4) / 1024 **2 ~ 356 Mb of float32 data
# and a rought test on the login node shows one slice taking about 9 minutes to load.

chunks = {"time": 744, "lat": 360, "lon": stripe_width}

# Read in the P Model time series variables
# - in each case, _open_ the dataset as an xarray MFDataset and then compute the subset
#   to be processed in the particular run

# ----------------------------------
# TEMPERATURE DATA
# ----------------------------------

# Find the temperature source files and open them as a dataset
temp_files = list((wfde_path / "Tair").rglob("*.nc"))
temp_source = xarray.open_mfdataset(temp_files, chunks=chunks)

# Create an indexing object as a tuple of slices that can be reused to get the variable
# subsets. Note that WFDE5 ends in 2019, hence the None.
idx_obj = (
    slice(
        np.where(temp_source.time.values == start_time)[0][0],
        None,
    ),
    slice(temp_source.sizes["lat"]),
    slice(lon_idx, lon_idx + stripe_width),
)

# Extract temperature data and convert to °C
temp_data = temp_source["Tair"][idx_obj].compute().data - 273.15

# Possible sources for atmospheric pressure
use_constant_patm = True

if use_constant_patm:
    # Use elevation derived atmospheric pressure
    # - Reduce to longitudinal slice and get patm
    elev_slice = asurf_data.data[tuple(list(idx_obj[1:]))]
    patm_slice = calc_patm(elev_slice)
    # - broadcast to shape of time series
    patm_data = np.broadcast_to(patm_slice, temp_data.shape)
else:
    # Find the atmospheric pressure source files and open them as a dataset
    patm_files = list((wfde_path / "PSurf").rglob("*.nc"))
    patm_source = xarray.open_mfdataset(patm_files, chunks=chunks)

    # Extract atmospheric pressure in Pa
    patm_data = patm_source["PSurf"][idx_obj].compute().data

# ----------------------------------
# PPFD DATA
# ----------------------------------

# Find the downwelling shortwave radiation files and open them as a dataset
swdown_files = list((wfde_path / "SWdown").rglob("*.nc"))
swdown_source = xarray.open_mfdataset(swdown_files, chunks=chunks)

# Extract shortwave downwelling radiation and convert to PPFD
ppfd_data = swdown_source["SWdown"][idx_obj].compute().data * 2.04

# ----------------------------------
# VPD DATA
# ----------------------------------

# Find the specific humidity files and open them as a dataset
qair_files = list((wfde_path / "Qair").rglob("*.nc"))
qair_source = xarray.open_mfdataset(qair_files, chunks=chunks)

# Extract specific humidity and convert to VPD: kg kg-1 to Pa.
qair_data = qair_source["SWdown"][idx_obj].compute().data
vpd_data = convert_sh_to_vpd(sh=qair_data, ta=temp_data, patm=patm_data)

# ----------------------------------
# FAPAR DATA
# ----------------------------------

# Load the daily fAPAR and then convert to subdaily
# NOTE - explicitly only getting files in C21 - problem with one file pre 2000
fapar_files = list(fapar_path.rglob("halfd_FPAR_20*.nc"))
fapar_source = xarray.open_mfdataset(fapar_files, chunks=chunks)

# Extract shortwave downwelling radiation and convert to PPFD
fapar_data = fapar_source["FPAR"][idx_obj].compute().data

# ----------------------------------
# CO2 DATA
# ----------------------------------

# Load the monthly co2 global time series and convert to a grid of subdaily values
co2_time_series = pd.read_csv(co2_path, comment="#")

# Get the datetime coordinates of the start of each month
month_start = (
    co2_time_series.year.astype(str)
    + "-"
    + co2_time_series.month.astype(str).str.zfill(2)
    + "-01"
)

# Convert to an array and resample to the hourly timesteps
hourly_co2 = (
    xarray.DataArray(
        co2_time_series["trend"], coords={"time": month_start.astype("datetime64[D]")}
    )
    .resample(time="1h")
    .ffill()
)

# Slice to the start time and endtime
hourly_co2 = hourly_co2[
    slice(
        np.where(hourly_co2.time.values == start_time)[0][0],
        np.where(hourly_co2.time.values == end_time)[0][0],
    )
]
