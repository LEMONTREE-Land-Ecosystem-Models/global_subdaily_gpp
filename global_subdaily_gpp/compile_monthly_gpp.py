"""This script takes the results files output by global_subdaily_models.py, which are
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
results_files = [p for (p, lon) in results_files]

start_time = time.time()
print("Starting script")

# Extract each month
# Something bizarre happens here with trying to open the full set of 720 files using
# open_mfdataset on the HPC. It very rapidly falls over with the message
# "AttributeError: NetCDF: Not a valid ID". Running a smaller set of files seems to
# succeed. For this reason, the loop below merges in a two step process.

results_blocks = np.split(np.array(results_files), 4)

for month in range(1, 13):

    def _preprocess(ds):
        return ds.sel(time=f"{year}-{month:02d}").drop_vars("local_time")

    blocks = []

    for idx, this_block in enumerate(results_blocks):
        block_results = xarray.open_mfdataset(
            this_block,
            parallel=True,
            preprocess=_preprocess,
            chunks={"lat": 360, "lon": 720, "time": 1},
            autoclose=True,
        )
        print(
            f"Year {year}, month {month:02d}, block {idx} opened"
            f"after {time.time() - start_time} seconds"
        )

        blocks.append(block_results.compute())

        print(
            f"Year {year}, month {month:02d}, block {idx} computed"
            f"after {time.time() - start_time} seconds"
        )

    results_month = xarray.merge(blocks)
    results_month.to_netcdf(months_dir / f"gpp_data_{year}_{month:02d}.nc")

    print(
        f"Year {year}, month {month:02d} written "
        f"after {time.time() - start_time} seconds"
    )
