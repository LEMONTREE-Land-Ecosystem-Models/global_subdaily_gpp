"""This script contains code to review the contents of an input data file created by
global_subdaily_gpp.py and visualise the data across time and latitudinal axes to check
data alignment."""

import matplotlib.pyplot as plt

import numpy as np
import pandas
import xarray
from pyrealm.pmodel import PModel, PModelEnvironment, FastSlowPModel, FastSlowScaler


data = xarray.open_dataset("results/inputs_data_163.nc")
global_pred = xarray.open_dataset("results/gpp_data_163.nc")
local_data = pandas.read_csv(
    "sites_comparison/giulia_data/dataOr_US-ARc#US-ARc_test_3a2005-2006.csv"
)

temp_data = data.temp
patm_data = data.patm
vpd_data = data.vpd
co2_data = data.co2
fapar_data = data.fapar
ppfd_data = data.ppfd


# Plot time series of one cell
fig, axes = plt.subplots(ncols=1, nrows=6)

for idx, dat in enumerate(
    [temp_data, patm_data, vpd_data, co2_data, fapar_data, ppfd_data]
):
    ax = axes[idx]
    ax.plot(dat[0 : (24 * 365), 284, 0])

plt.tight_layout()
plt.show()


# Plot latitudinal profile of data availablility
fig, axes = plt.subplots(ncols=1, nrows=6)

for idx, dat in enumerate(
    [temp_data, patm_data, vpd_data, co2_data, fapar_data, ppfd_data]
):
    ax = axes[idx]
    ax.plot(dat.count(dim=["time", "lon"]))

plt.tight_layout()
plt.show()


# Adjust the UTC time to be the local time
data = data.assign_coords(
    time=data.time.data + np.timedelta64(round(((-98.25 / 180) * 12) * 60), "m")
)

# Reduce to focal cell and time slice of local data
local_datetime = pandas.to_datetime(local_data["TIME"])

cell_data = data.sel(
    lat=35.75,
    lon=-98.25,
    time=slice(local_datetime.min(), local_datetime.max()),
)


# Extract variables and forward fill to remove NAs
cell_temp = cell_data.temp
cell_patm = cell_data.patm
cell_vpd = cell_data.vpd
cell_co2 = cell_data.co2
cell_fapar = cell_data.fapar
cell_ppfd = cell_data.ppfd

# Get the P Model environment
cell_pm_env = PModelEnvironment(
    tc=cell_temp.data, patm=cell_patm.data, vpd=cell_vpd.data, co2=cell_co2.data
)

cell_standard_pmod = PModel(cell_pm_env, kphio=1 / 8)

cell_standard_pmod.estimate_productivity(fapar=cell_fapar.data, ppfd=cell_ppfd.data)

fsscalar = FastSlowScaler(cell_data.time.data)
fsscalar.set_window(
    window_center=np.timedelta64(12, "h"),
    half_width=np.timedelta64(1, "h"),
)

cell_subdaily_pmod = FastSlowPModel(
    env=cell_pm_env,
    fs_scaler=fsscalar,
    ppfd=cell_ppfd.data,
    fapar=cell_fapar.data,
)


# Extract hourly predictions from half hourly local data
local_hourly = local_data.loc[local_datetime.dt.minute == 0]

# Plot global against local
plt.plot(local_hourly["GPPp"].values * 12.011)
plt.plot(cell_standard_pmod.gpp)
plt.show()

plt.plot(local_hourly["GPPpOpt"].values * 12.011)
plt.plot(cell_subdaily_pmod.gpp)
plt.show()

fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True)
ax1.plot(cell_data.time.data, cell_standard_pmod.gpp)
ax2.plot(cell_data.time.data, cell_subdaily_pmod.gpp)
plt.show()


# Compare global values calculated here against global from HPC

# Adjust the UTC time to be the local time
global_pred = global_pred.assign_coords(
    time=global_pred.time.data + np.timedelta64(round(((-98.25 / 180) * 12) * 60), "m")
)

global_pred_cell_data = global_pred.sel(
    lat=35.75,
    lon=-98.25,
    time=slice(local_datetime.min(), local_datetime.max()),
)

# Yeah something wrong with the global implementation.
n = 300
plt.plot(cell_data.time.data[0:n], cell_standard_pmod.gpp[0:n])
plt.plot(cell_data.time.data[0:n], cell_subdaily_pmod.gpp[0:n])
plt.plot(global_pred_cell_data.time.data[0:n], global_pred_cell_data.subdaily_gpp[0:n])

plt.show()


# Crop time offset data to complete days
global_trim_data = data.sel(time=slice("2005-01-01", "2007-01-01"))

# Extract variables and forward fill to remove NAs
gt_cell_temp = global_trim_data.temp
gt_cell_patm = global_trim_data.patm
gt_cell_vpd = global_trim_data.vpd
gt_cell_co2 = global_trim_data.co2
gt_cell_fapar = global_trim_data.fapar
gt_cell_ppfd = global_trim_data.ppfd

# Get the P Model environment
global_pm_env = PModelEnvironment(
    tc=gt_cell_temp.data,
    patm=gt_cell_patm.data,
    vpd=gt_cell_vpd.data,
    co2=gt_cell_co2.data,
)

global_standard_pmod = PModel(global_pm_env, kphio=1 / 8)

global_standard_pmod.estimate_productivity(
    fapar=gt_cell_fapar.data, ppfd=gt_cell_ppfd.data
)

global_fsscalar = FastSlowScaler(global_trim_data.time.data)
global_fsscalar.set_window(
    window_center=np.timedelta64(12, "h"),
    half_width=np.timedelta64(1, "h"),
)

global_subdaily_pmod = FastSlowPModel(
    env=global_pm_env,
    fs_scaler=fsscalar,
    ppfd=gt_cell_ppfd.data,
    fapar=gt_cell_fapar.data,
)


# Yeah something wrong with the global implementation.
n = 300
plt.plot(
    cell_data.time.data[0:n],
    cell_standard_pmod.gpp[0:n],
    linewidth=3,
)
plt.plot(
    global_trim_data.time.data[0:n],
    global_standard_pmod.gpp[0:n, 251, 0],
)

plt.plot(
    cell_data.time.data[0:n],
    cell_subdaily_pmod.gpp[0:n],
    linewidth=3,
)

plt.plot(
    global_trim_data.time.data[0:n],
    global_subdaily_pmod.gpp[0:n, 251, 0],
)

plt.plot(
    global_pred_cell_data.time.data[0:n],
    global_pred_cell_data.subdaily_gpp[0:n],
)


plt.show()
