"""Fitting the subdaily P Model to global data.

This script is to be used as part of an array job on the HPC to iterate over blocks of
cells from WFDE5 data in order to fit the subdaily P Model at global scale and 0.5°
resolution.

The PBS_ARRAY_INDEX environment variable changes within each subjob and is used to
control which longitudinal slice is analysed by each subjob.
"""

import os
from pathlib import Path
import time

import xarray
import numpy as np
import pandas as pd

from pyrealm.hygro import convert_sh_to_vpd
from pyrealm.pmodel.functions import calc_patm
from pyrealm.pmodel import PModel, PModelEnvironment, FastSlowPModel, FastSlowScaler


# ----------------------------------
# DATA PATHS
# ----------------------------------
script_start = time.time()

# WFDE5 provides TEMP, PATM, VPD, PPFD
wfde_path = Path("/rds/general/project/lemontree/live/source/wfde5/wfde5_v2")
# SNU FAPAR downsampled from 0.05° to 0.5° gives FAPAR
fapar_path = Path("/rds/general/project/lemontree/live/derived/SNU_Ryu/FPAR2_05d")
# NOAA global sea level time series provides CO2
co2_path = Path("/rds/general/project/lemontree/live/source/NOAA_CO2/co2_mm_gl.csv")
# Output directory
output_path = Path(
    "/rds/general/project/lab-prentice-realm-data/live/splash_CRU/"
    "global_subdaily_gpp/results"
)


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
    slice(asurf_data.sizes["lat"]),
    slice(lon_idx, lon_idx + stripe_width),
)

# Retain the temperature xarray to use the coordinates for outputs.
temp_xarray = temp_source["Tair"][idx_obj]

# Extract temperature data and convert to °C
temp_data = temp_xarray.compute().data - 273.15

# Remove very cold cells
temp_data[temp_data < -25.0] = np.nan

# ----------------------------------
# PATM DATA
# ----------------------------------

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
# VPD DATA
# ----------------------------------

# Find the specific humidity files and open them as a dataset
qair_files = list((wfde_path / "Qair").rglob("*.nc"))
qair_source = xarray.open_mfdataset(qair_files, chunks=chunks)

# Extract specific humidity and convert to VPD: kg kg-1 to Pa.
qair_data = qair_source["Qair"][idx_obj].compute().data
# Function takes pressure in kPa and returns kPa
vpd_data = convert_sh_to_vpd(sh=qair_data, ta=temp_data, patm=patm_data / 1000) * 1000

# Set negative values to zero
vpd_data = np.clip(vpd_data, 0, np.inf)

# ----------------------------------
# PPFD DATA
# ----------------------------------

# Find the downwelling shortwave radiation files and open them as a dataset
swdown_files = list((wfde_path / "SWdown").rglob("*.nc"))
swdown_source = xarray.open_mfdataset(swdown_files, chunks=chunks)

# SWDown has a different start datetime to other WFDE5 variables, so needs different
# indexing
idx_obj = (
    slice(
        np.where(swdown_source.time.values == start_time)[0][0],
        None,
    ),
    slice(asurf_data.sizes["lat"]),
    slice(lon_idx, lon_idx + stripe_width),
)

# Extract shortwave downwelling radiation and convert to PPFD
ppfd_data = swdown_source["SWdown"][idx_obj].compute().data * 2.04

# ----------------------------------
# FAPAR DATA
# ----------------------------------

# Load the daily fAPAR and then convert to subdaily
fapar_files = list(fapar_path.rglob("*.nc"))
fapar_source = xarray.open_mfdataset(fapar_files, chunks=chunks)

# Forward fill FAPAR to hourly sampling
fapar_hourly = fapar_source["FPAR"].resample(time="1h").ffill()

# Get the indices on the hourly resample
idx_obj_fapar = (
    slice(
        np.where(fapar_hourly.time.values == start_time)[0][0],
        np.where(fapar_hourly.time.values == end_time)[0][0],
    ),
    slice(asurf_data.sizes["lat"]),
    slice(lon_idx, lon_idx + stripe_width),
)

# Extract the data
fapar_data = fapar_hourly[idx_obj_fapar].compute().data

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

# Broadcast the data to the correct shape
co2_data = np.broadcast_to(hourly_co2.data[:, np.newaxis, np.newaxis], temp_data.shape)

print(
    f"""
Data loading finished after {time.time() - script_start} seconds:
- lon_idx = {lon_idx}
- temp_data.shape = {temp_data.shape}
- patm_data.shape = {patm_data.shape}
- vpd_data.shape = {vpd_data.shape}
- co2_data.shape = {co2_data.shape}
- fapar_data.shape = {fapar_data.shape}
- ppfd_data.shape = {ppfd_data.shape}
"""
)

# ----------------------------------
# MODEL FITTING
# ----------------------------------

# Loop over the loaded longitudinal bands, correcting the datetimes from UTC to local
# time for calculating subdaily representative times

utc_times = temp_xarray.time.values
lon_bands = np.arange(lon_idx, lon_idx + stripe_width)


results = []

for data_idx, this_lon_idx in enumerate(lon_bands):
    # Indexer for the longitudinal band
    lidx = (
        slice(temp_xarray.sizes["time"]),
        slice(asurf_data.sizes["lat"]),
        slice(data_idx, data_idx + 1),
    )

    # Xarray coordinates for the band
    band_coords = temp_xarray[lidx].coords

    # Get the P Model environment
    pm_env = PModelEnvironment(
        tc=temp_data[lidx], patm=patm_data[lidx], vpd=vpd_data[lidx], co2=co2_data[lidx]
    )

    # Print out a data summary for the photosynthetic environment
    pm_env.summarize()

    # Fit the standard P Model
    standard_pmod = PModel(pm_env, kphio=1 / 8)
    standard_pmod.estimate_productivity(fapar=fapar_data[lidx], ppfd=ppfd_data[lidx])

    # Print out a summary for the standard model
    standard_pmod.summarize()

    # Set up the local times - the lon idx in 0-719 defines 0.5° longitudinal bands,
    # each of width is 2 minutes wide (24 / (720) * 60 = 2.0). So, find the local time
    # for the left hand edge using this_lon_idx and add 1 minute to get the band centre.

    local_time_delta = (24 * (this_lon_idx / 720) - 12) * 60 + 1
    local_time_delta = np.timedelta64(int(local_time_delta), "m")
    local_time = utc_times + local_time_delta

    # Set a half hourly window around noon - with hourly data this is actually just
    # picking the noon value.
    fsscaler = FastSlowScaler(local_time)
    fsscaler.set_window(
        window_center=np.timedelta64(12, "h"),
        half_width=np.timedelta64(1, "h"),
    )

    # Fit the subdaily P Model
    fs_pmod = FastSlowPModel(
        env=pm_env,
        fs_scaler=fsscaler,
        fapar=fapar_data[lidx],
        ppfd=ppfd_data[lidx],
        alpha=1 / 15,
        kphio=1 / 8,
    )

    # Format and store the GPP data
    res = xarray.Dataset(
        {
            "standard_gpp": xarray.DataArray(temp_data, coords=band_coords),
            "subdaily_gpp": xarray.DataArray(temp_data, coords=band_coords),
        }
    )
    results.append(res)

    print(
        f"Models fitted on lon_idx {this_lon_idx} after "
        f"{time.time() - script_start} seconds"
    )

# Concatenate the results along the longitude axis and export
out_data = xarray.concat(results, dim="lon")
out_data.to_netcdf(output_path / f"gpp_data_{array_index}.nc")

print(f"Script completed after {time.time() - script_start} seconds")
