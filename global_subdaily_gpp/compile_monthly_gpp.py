"""This script takes the results files output by global_subdaily_models.py, which are
0.5° longitudinal slices of a 20 year time series, and converts back to regular lat long
grids. It partitions the time series into more manageable monthly output files.

It is written to be run as an array job, with each subjob handling a subset of the time
series.
"""

import os
from pathlib import Path

import xarray


results_dir = Path(
    "/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"
    "global_subdaily_gpp/results"
)

months_dir = Path(
    "/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"
    "global_subdaily_gpp/monthly_gpp"
)

# ----------------------------------
# IDENTIFY TIME BLOCK
# ----------------------------------
array_index = int(os.environ["PBS_ARRAY_INDEX"])
year = array_index + 2000

# ----------------------------------
# BUILD A MULTIPLE FILE DATASET OF THE 0.5° BAND RESULTS
# ----------------------------------

# Find the 0.5° band output files
results_files = list(results_dir.rglob("*.nc"))


# Extract each month

for month in range(1, 13):

    def _preprocess(ds):
        return ds.sel(time=f"{year}-{month:02d}").drop_vars("local_time")

    results_month = xarray.open_mfdataset(
        results_files,
        parallel=True,
        preprocess=_preprocess,
        chunks={"lat": 360, "lon": 720, "time": 1},
    )

    results_month.to_netcdf(months_dir / f"gpp_data_{year}_{month:02d}.nc")
