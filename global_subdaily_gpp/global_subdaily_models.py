"""Fitting the subdaily P Model to global data.

This script is to be used as part of an array job on the HPC to iterate over blocks of
cells from WFDE5 data in order to fit the subdaily P Model at global scale and 0.5°
resolution.

The PBS_ARRAY_INDEX environment variable changes within each subjob and is used to
control which block of longitudinal slices is tackled by the subjob. The N_LON_SLICES
environment variable is used to set how many blocks the 720 longitudinal bands are
divided into to be handled.

The timing and data requirements are roughly that loading and calculating the outputs
for one longitudinal band (360 0.5° lat cells x 1 0.5° long cell x ~ 175200 hours of
data over 20 years) runs at about 15 minites data loading and aligning and 5 minutes
modelling with peak memory usage of 55GB. It is probably more efficient to load more
data in one go, but multiples of 55GB don't really scale on the general/throughput
nodes, so this also iterates the data loading.

So: using N_LON_SLICES=30 gives 30 subjobs, each tackling 24 x 0.5° longitudinal bands
and each running for around 24 * 20 / 60 = 8 hours, with some extra for variable
runspeeds.

By default, the number of longitudinal slices is set to the -179.75 to 179.75 range of
the input grids. However errors occurring on the HPC can leave blocks incomplete. To
make it easier to fill these gaps, a file of target longitude bands can be provided
using the TARGET_LON_BANDS environment variable. This file should simply contain the
longitude values of the bands to run, one per line.

If the WRITE_PMODEL_INPUTS environment variable is set to any value, the assembled,
aligned and filled dataset for longitudinal slices will be written out to allow
debugging, but these are not small!
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
end_time = np.datetime64("2019-12-31 23:59")

# ----------------------------------
# ELEVATION DATA
# ----------------------------------

# Read in the elevation data
asurf_data = xarray.load_dataarray(wfde_path / "Elev" / "ASurf_WFDE5_CRU_v2.0.nc")

# ----------------------------------
# IDENTIFY LONGITUDINAL BLOCKS
# ----------------------------------

# Use longitudinal stripes to control for timing: data is sampled at UTC, so local time
# varies with longitude and needs to be corrected. The blocksize sets the number of
# those vertical stripes to be loaded

# The array index gives the chunk number given the number of chunks, defaulting to a
# single longitudinal band per chunk
n_chunks = int(os.environ.get("N_LON_SLICES", 720))

target_lon_bands = os.environ.get("TARGET_LON_BANDS")
if target_lon_bands is None:
    # Divide the longitudinal coordinates of the input grid into chunks
    lon_chunks = np.array_split(asurf_data.lon.data, n_chunks)
else:
    # Divide the actual longitudinal values provided in an input file into chunks
    lon_chunks = np.array_split(np.loadtxt(target_lon_bands), n_chunks)


# PBS job array indices start at 1
array_index = int(os.environ["PBS_ARRAY_INDEX"])
lon_vals = lon_chunks[array_index]

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

chunks = {"time": 744, "lat": 360, "lon": 1}

# Read in the P Model time series variables
# - in each case, _open_ the dataset as an xarray MFDataset and then compute the subset
#   to be processed in the particular run

# ----------------------------------
# BUILD THE MULTIPLE FILE DATASET LINKS TO THE SOURCE FILES TO BE RE-USED ACROSS BANDS
# ----------------------------------

# Find the hourly WDFE temperature source files and open them as a dataset
temp_files = list((wfde_path / "Tair").rglob("*.nc"))
temp_source = xarray.open_mfdataset(temp_files, chunks=chunks)

# Find the hourly WDFE atmospheric pressure source files and open them as a dataset
patm_files = list((wfde_path / "PSurf").rglob("*.nc"))
patm_source = xarray.open_mfdataset(patm_files, chunks=chunks)

# Find the hourly WDFE specific humidity files and open them as a dataset
qair_files = list((wfde_path / "Qair").rglob("*.nc"))
qair_source = xarray.open_mfdataset(qair_files, chunks=chunks)

# Find the hourly WDFE downwelling shortwave radiation files and open them as a dataset
swdown_files = list((wfde_path / "SWdown").rglob("*.nc"))
swdown_source = xarray.open_mfdataset(swdown_files, chunks=chunks)

# Find the daily SNU FAPAR files and open them as a dataset
fapar_files = list(fapar_path.rglob("*.nc"))
fapar_source = xarray.open_mfdataset(fapar_files, chunks=chunks)

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
        co2_time_series["trend"],
        coords={"time": month_start.astype("datetime64[D]")},
    )
    .resample(time="1h")
    .ffill()
)

# Slice to the start time and endtime
hourly_co2 = hourly_co2.sel({"time": slice(start_time, end_time)})

# ----------------------------------
# LOOP OVER THE LONGITUDINAL BANDS
# ----------------------------------

for this_lon_val in lon_vals:
    # Create an indexing dictionary to subset the time and longitude axes of the WFDE
    # datasets - note the use of a list in 'lon' to _retain_ that as a singleton axis.
    wfde_slices = {
        "time": slice(start_time, end_time),
        "lon": [this_lon_val],
    }

    # ----------------------------------
    # TEMPERATURE DATA
    # ----------------------------------

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
        elev_slice = asurf_data.sel({"lon": [this_lon_val]})
        patm_slice = calc_patm(elev_slice.data)
        # - broadcast to shape of time series
        patm_data = xarray.DataArray(
            np.broadcast_to(patm_slice, temp_data.shape), coords=temp_data.coords
        )
    else:
        # Extract atmospheric pressure in Pa
        patm_data = patm_source["PSurf"].sel(wfde_slices).compute()

    # ----------------------------------
    # VPD DATA
    # ----------------------------------

    # Extract specific humidity and convert to VPD: kg kg-1 to Pa.
    qair_data = qair_source["Qair"].sel(wfde_slices).compute().data

    # Function takes pressure in kPa and returns kPa
    vpd_data = (
        convert_sh_to_vpd(sh=qair_data, ta=temp_data, patm=patm_data / 1000) * 1000
    )

    # Set negative values to zero
    vpd_data = np.clip(vpd_data, 0, np.inf)

    # Convert to xarray
    vpd_data = xarray.DataArray(vpd_data, coords=temp_data.coords)

    # ----------------------------------
    # PPFD DATA
    # ----------------------------------

    # Extract shortwave downwelling radiation and convert to PPFD
    ppfd_data = swdown_source["SWdown"].sel(wfde_slices).compute() * 2.04

    # ----------------------------------
    # FAPAR DATA
    # ----------------------------------

    # Forward fill FAPAR to hourly sampling
    fapar_hourly = fapar_source["FPAR"].resample(time="1h").ffill()

    # Extract the data
    fapar_data = fapar_hourly.sel(wfde_slices).compute()

    # ----------------------------------
    # CO2 DATA
    # ----------------------------------

    # Broadcast the data to the correct shape and convert to xarray
    co2_data = np.broadcast_to(
        hourly_co2.data[:, np.newaxis, np.newaxis], temp_data.shape
    )

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

    # ----------------------------------
    # COMPILE INTO A SINGLE DATASET AND FORWARD FILL MISSING DATA
    # ----------------------------------

    this_lon_inputs = xarray.Dataset(
        {
            "temp": temp_data,
            "patm": patm_data,
            "vpd": vpd_data,
            "co2": co2_data,
            "fapar": fapar_data,
            "ppfd": ppfd_data,
        }
    )

    this_lon_inputs = this_lon_inputs.ffill(dim="time")

    print(f"Dataset built and gaps filled after {time.time() - script_start} seconds")

    # ----------------------------------
    # SAVE MODEL INPUTS IF REQUESTED
    # ----------------------------------

    if os.environ.get("WRITE_PMODEL_INPUTS"):
        out = this_lon_inputs.to_netcdf(output_path / f"inputs_data_{this_lon_val}.nc")

        print(f"Input data saved after {time.time() - script_start} seconds")

    # ----------------------------------
    # MODEL FITTING
    # ----------------------------------

    # Correcting the datetimes from UTC to local time for calculating subdaily
    # representative times and run the models.

    # Set up the local times - the lon idx in 0-719 defines 0.5° longitudinal bands,
    # each of width is 2 minutes wide (24 / (720) * 60 = 2.0). So, find the local time
    # for the left hand edge using this_lon_idx and add 1 minute to get the band centre.
    local_time_delta = ((this_lon_val / 180) * 12) * 60 + 1
    local_time_delta = np.timedelta64(round(local_time_delta), "m")
    local_time = this_lon_inputs.time + local_time_delta

    # Assign the correct local time coordinates for this longitude slice
    this_lon_inputs = this_lon_inputs.assign_coords(time=local_time)

    # Moving to local times shifts the datetimes from complete days to partial days,
    # which are not currently supported by the FastSlowScaler. So need to reduce the
    # analysis to complete days by trimming a day from each end to remove partials.
    this_lon_inputs = this_lon_inputs.sel(
        time=slice(
            start_time + np.timedelta64(1, "D"),
            end_time - np.timedelta64(1, "D"),
        )
    )

    # Get the P Model environment
    pm_env = PModelEnvironment(
        tc=this_lon_inputs["temp"].data,
        patm=this_lon_inputs["patm"].data,
        vpd=this_lon_inputs["vpd"].data,
        co2=this_lon_inputs["co2"].data,
    )

    # Print out a data summary for the photosynthetic environment
    pm_env.summarize()

    # Fit the standard P Model
    standard_pmod = PModel(pm_env, kphio=1 / 8)
    standard_pmod.estimate_productivity(
        fapar=this_lon_inputs["fapar"].data,
        ppfd=this_lon_inputs["ppfd"].data,
    )

    # Print out a summary for the standard model
    standard_pmod.summarize()

    # Set a half hourly window around noon - with hourly data this is actually just
    # picking the noon value.
    fsscaler = FastSlowScaler(this_lon_inputs.time.data)
    fsscaler.set_window(
        window_center=np.timedelta64(12, "h"),
        half_width=np.timedelta64(1, "h"),
    )

    # Fit the subdaily P Model
    subdaily_pmod = FastSlowPModel(
        env=pm_env,
        fs_scaler=fsscaler,
        handle_nan=True,
        fapar=this_lon_inputs["fapar"].data,
        ppfd=this_lon_inputs["ppfd"].data,
        alpha=1 / 15,
        kphio=1 / 8,
    )

    print(
        f"Models fitted on lon {this_lon_val} after "
        f"{time.time() - script_start} seconds"
    )

    # Format and store the GPP data
    res = xarray.Dataset(
        {
            "standard_gpp": xarray.DataArray(
                standard_pmod.gpp,
                dims=["time", "lat", "lon"],
                coords=this_lon_inputs.coords,
            ).astype(np.float32),
            "subdaily_gpp": xarray.DataArray(
                subdaily_pmod.gpp,
                dims=["time", "lat", "lon"],
                coords=this_lon_inputs.coords,
            ).astype(np.float32),
        }
    )

    print(
        f"Data added to results on lon {this_lon_val} after "
        f"{time.time() - script_start} seconds"
    )

    res.to_netcdf(output_path / f"gpp_data_{this_lon_val}.nc")

    print(f"Data saved and script completed after {time.time() - script_start} seconds")
