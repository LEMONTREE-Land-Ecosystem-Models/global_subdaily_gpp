"""
This script takes the results files output by global_subdaily_models.py, which are
0.5° longitudinal slices of a 20 year time series, and converts back to regular lat long
grids. It partitions the time series into more manageable monthly output files.

It is written to be run as an array job, with each subjob handling a subset of the time
series. 
"""

import os
from pathlib import Path
import time

import xarray
import numpy as np

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

# Find the 0.5° band output files and sort by longitude
results_files = list(results_dir.rglob("gpp_data*.nc"))
results_files = [(p, float(p.name[9:-3])) for p in results_files]
results_files.sort(key=lambda x: x[1])


start_time = time.time()
print("Starting script")

# Extract each month

for month in range(1, 13):
    data = []

    for this_file, lon in results_files:
        # Open the source file using a with block to close cleanly and append the data
        # for the month, dropping the local time variable, which changes across each
        # longitudinal band.
        with xarray.open_dataset(this_file) as ds:
            data.append(ds.sel(time=f"{year}-{month:02d}").drop_vars("local_time"))

        print(
            f"Year {year}, month {month:02d}, lon {lon} loaded"
            f"after {time.time() - start_time} seconds"
        )

    results_month = xarray.merge(data)
    results_month.to_netcdf(months_dir / f"gpp_data_{year}_{month:02d}.nc")

    print(
        f"Year {year}, month {month:02d} written "
        f"after {time.time() - start_time} seconds"
    )
