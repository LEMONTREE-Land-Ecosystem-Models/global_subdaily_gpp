"""This script is used to review the contents of an input data file and visualise the
data across time and latitudinal axes to check data alignment."""

import matplotlib.pyplot as plt

import numpy as np
import xarray
from pyrealm.pmodel import PModel, PModelEnvironment, FastSlowPModel, FastSlowScaler


data = xarray.open_dataset("global_subdaily_gpp/results/inputs_data_408.nc")

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


# Reduce to focal cell
cell_data = data.sel(lat=61.75, lon=24.25)

# Extract variables and forward fill to remove NAs
cell_temp = cell_data.temp.ffill(dim="time")
cell_patm = cell_data.patm.ffill(dim="time")
cell_vpd = cell_data.vpd.ffill(dim="time")
cell_co2 = cell_data.co2.ffill(dim="time")
cell_fapar = cell_data.fapar.ffill(dim="time")
cell_ppfd = cell_data.ppfd.ffill(dim="time")


# Get the P Model environment
pm_env = PModelEnvironment(
    tc=cell_temp.data, patm=cell_patm.data, vpd=cell_vpd.data, co2=cell_co2.data
)

standard_pmod = PModel(pm_env, kphio=1 / 8)

standard_pmod.estimate_productivity(fapar=cell_fapar.data, ppfd=cell_ppfd.data)


fsscalar = FastSlowScaler(cell_data.time.data)
fsscalar.set_window(
    window_center=np.timedelta64(12, "h"),
    half_width=np.timedelta64(1, "h"),
)


subdaily_pmod = FastSlowPModel(
    env=pm_env,
    fs_scaler=fsscalar,
    ppfd=cell_ppfd.data,
    fapar=cell_fapar.data,
)

fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True)
ax1.plot(cell_data.time.data, standard_pmod.gpp)
ax2.plot(cell_data.time.data, subdaily_pmod.gpp)
plt.show()
