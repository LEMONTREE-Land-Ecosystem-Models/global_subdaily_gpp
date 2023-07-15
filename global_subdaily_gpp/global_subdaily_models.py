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
fapar_path = Path("/rds/general/project/lemontree/live/derived/SNU_Ryu/FPAR_05d")
# NOAA global sea level time series provides CO2
co2_path = Path("/rds/general/project/lemontree/live/source/NOAA_CO2/co2_mm_gl.csv")
# Output directory
output_path = Path(
    "/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"
    "global_subdaily_gpp/results"
)


start_time = np.datetime64("2000-01-01 00:00")
end_time = np.datetime64("2000-12-31 23:59")

# ----------------------------------
# ELEVATION DATA
# ----------------------------------

# Read in the elevation data
asurf_data = xarray.load_dataarray(wfde_path / "Elev" / "ASurf_WFDE5_CRU_v2.0.nc")

# Use longitudinal stripes to control for timing: data is sampled at UTC, so local time
# varies with longitude and needs to be corrected. The blocksize sets the number of
# those vertical stripes to be loaded
array_index = int(os.environ["PBS_ARRAY_INDEX"])
stripe_width = int(os.environ.get("N_LON_SLICES", 1))
lon_vals = asurf_data.lon[array_index : (array_index + stripe_width)].data

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

# Create an indexing dictionary to subset the time and longitude axes of the WFDE
# datasets - note the use of a list in 'lon' to _retain_ that as a singleton axis.
wfde_slices = {
    "time": slice(start_time, end_time),
    "lon": lon_vals,
}


# Extract temperature data and convert to °C
temp_data = temp_source["Tair"].sel(wfde_slices).compute() - 273.15

# Remove very cold cells
temp_data = temp_data.where(temp_data >= -25.0, np.nan)

# ----------------------------------
# PATM DATA
# ----------------------------------

# Possible sources for atmospheric pressure
use_constant_patm = True

if use_constant_patm:
    # Use elevation derived atmospheric pressure
    # - Reduce to longitudinal slice and get patm
    elev_slice = asurf_data.sel({"lon": lon_vals})
    patm_slice = calc_patm(elev_slice.data)
    # - broadcast to shape of time series
    patm_data = xarray.DataArray(
        np.broadcast_to(patm_slice, temp_data.shape), coords=temp_data.coords
    )
else:
    # Find the atmospheric pressure source files and open them as a dataset
    patm_files = list((wfde_path / "PSurf").rglob("*.nc"))
    patm_source = xarray.open_mfdataset(patm_files, chunks=chunks)

    # Extract atmospheric pressure in Pa
    patm_data = patm_source["PSurf"].sel(wfde_slices).compute()

# ----------------------------------
# VPD DATA
# ----------------------------------

# Find the specific humidity files and open them as a dataset
qair_files = list((wfde_path / "Qair").rglob("*.nc"))
qair_source = xarray.open_mfdataset(qair_files, chunks=chunks)

# Extract specific humidity and convert to VPD: kg kg-1 to Pa.
qair_data = qair_source["Qair"].sel(wfde_slices).compute().data
# Function takes pressure in kPa and returns kPa
vpd_data = convert_sh_to_vpd(sh=qair_data, ta=temp_data, patm=patm_data / 1000) * 1000

# Set negative values to zero
vpd_data = np.clip(vpd_data, 0, np.inf)

# Convert to xarray
vpd_data = xarray.DataArray(vpd_data, coords=temp_data.coords)

# ----------------------------------
# PPFD DATA
# ----------------------------------

# Find the downwelling shortwave radiation files and open them as a dataset
swdown_files = list((wfde_path / "SWdown").rglob("*.nc"))
swdown_source = xarray.open_mfdataset(swdown_files, chunks=chunks)


# Extract shortwave downwelling radiation and convert to PPFD
ppfd_data = swdown_source["SWdown"].sel(wfde_slices).compute() * 2.04

# ----------------------------------
# FAPAR DATA
# ----------------------------------

# Load the daily fAPAR and then convert to subdaily
fapar_files = list(fapar_path.rglob("*.nc"))
fapar_source = xarray.open_mfdataset(fapar_files, chunks=chunks)

# Forward fill FAPAR to hourly sampling
fapar_hourly = fapar_source["FPAR"].resample(time="1h").ffill()

# Extract the data
fapar_hourly = fapar_hourly.rename({"longitude": "lon", "latitude": "lat"})
fapar_data = fapar_hourly.sel(wfde_slices).compute()

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
hourly_co2 = hourly_co2.sel({"time": slice(start_time, end_time)})

# Broadcast the data to the correct shape and convert to xarray
co2_data = np.broadcast_to(hourly_co2.data[:, np.newaxis, np.newaxis], temp_data.shape)

co2_data = xarray.DataArray(co2_data, coords=temp_data.coords)

print(
    f"""
Data loading finished after {time.time() - script_start} seconds:
- lon_idx = {lon_vals}
- temp_data.shape = {temp_data.shape}
- patm_data.shape = {patm_data.shape}
- vpd_data.shape = {vpd_data.shape}
- co2_data.shape = {co2_data.shape}
- fapar_data.shape = {fapar_data.shape}
- ppfd_data.shape = {ppfd_data.shape}
"""
)

if os.environ.get("WRITE_PMODEL_INPUTS"):
    out = xarray.Dataset(
        {
            "temp": temp_data,
            "patm": patm_data,
            "vpd": vpd_data,
            "co2": co2_data,
            "fapar": fapar_data,
            "ppfd": ppfd_data,
        }
    ).to_netcdf(output_path / f"inputs_data_{array_index}.nc")

# ----------------------------------
# MODEL FITTING
# ----------------------------------

# Loop over the loaded longitudinal bands, correcting the datetimes from UTC to local
# time for calculating subdaily representative times

utc_times = temp_data.time.values

results = []

for this_lon in lon_vals:
    # Get an indexing dictionary for the longitudinal band
    lon_sel = {"lon": this_lon}

    # Xarray coordinates for the band
    band_coords = temp_data.sel(lon_sel).coords

    # Get the P Model environment
    pm_env = PModelEnvironment(
        tc=temp_data.data.sel(lon_sel),
        patm=patm_data.data.sel(lon_sel),
        vpd=vpd_data.data.sel(lon_sel),
        co2=co2_data.data.sel(lon_sel),
    )

    # Print out a data summary for the photosynthetic environment
    pm_env.summarize()

    # Fit the standard P Model
    standard_pmod = PModel(pm_env, kphio=1 / 8)
    standard_pmod.estimate_productivity(
        fapar=fapar_data.data.sel(lon_sel),
        ppfd=ppfd_data.data.sel(lon_sel),
    )

    # Print out a summary for the standard model
    standard_pmod.summarize()

    # Set up the local times - the lon idx in 0-719 defines 0.5° longitudinal bands,
    # each of width is 2 minutes wide (24 / (720) * 60 = 2.0). So, find the local time
    # for the left hand edge using this_lon_idx and add 1 minute to get the band centre.

    local_time_delta = ((this_lon / 180) * 12) * 60
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
    subdaily_pmod = FastSlowPModel(
        env=pm_env,
        fs_scaler=fsscaler,
        fapar=fapar_data.data.sel(lon_sel),
        ppfd=ppfd_data.data(lon_sel),
        alpha=1 / 15,
        kphio=1 / 8,
    )

    # Format and store the GPP data
    res = xarray.Dataset(
        {
            "standard_gpp": xarray.DataArray(standard_pmod.gpp, coords=band_coords),
            "subdaily_gpp": xarray.DataArray(subdaily_pmod.gpp, coords=band_coords),
        }
    )
    results.append(res)

    print(
        f"Models fitted on lon {this_lon} after "
        f"{time.time() - script_start} seconds"
    )

# Concatenate the results along the longitude axis and export
out_data = xarray.concat(results, dim="lon")
out_data.to_netcdf(output_path / f"gpp_data_{array_index}.nc")

print(f"Script completed after {time.time() - script_start} seconds")
