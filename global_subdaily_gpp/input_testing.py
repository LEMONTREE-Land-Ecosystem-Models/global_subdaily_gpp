"""This script is used to review the contents of an input data file and visualise the
data across time and latitudinal axes to check data alignment."""

import matplotlib.pyplot as plt

import xarray
from pyrealm.pmodel import PModel, PModelEnvironment, FastSlowPModel, FastSlowScaler


data = xarray.open_dataset("global_subdaily_gpp/results/inputs_data_360.nc")

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


# Get the P Model environment
pm_env = PModelEnvironment(
    tc=temp_data.data,
    patm=patm_data.data,
    vpd=vpd_data.data,
    co2=co2_data.data,
)

standard_pmod = PModel(pm_env, kphio=1 / 8)


standard_pmod.estimate_productivity(fapar=fapar_data.data, ppfd=ppfd_data.data)
