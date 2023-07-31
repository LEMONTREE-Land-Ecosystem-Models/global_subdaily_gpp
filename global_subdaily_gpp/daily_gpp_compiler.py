from pathlib import Path

import xarray
import numpy as np

results_dir = Path(
    "/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"
    "global_subdaily_gpp/results"
)

soil_beta_path = Path(
    "/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"
    "global_splash_cru/soil_moisture_grids"
)


# Find the 0.5° band output files and sort by longitude
results_files = list(results_dir.rglob("gpp_data*.nc"))
results_files = [(p, float(p.name[9:-3])) for p in results_files]
results_files.sort(key=lambda x: x[1])
results_files = [p for (p, lon) in results_files]

# Stores for processed data
standard_daily = []
subdaily_daily = []
subdaily_with_beta_daily = []

# Open soil beta penalty dataset from yearly files
soil_beta = xarray.open_mfdataset(soil_beta_path.glob("soil_moisture*.nc"))

for each_file in results_files[350:370]:
    # Loop over the longitudinal files, appending daily summaries for each band to the
    # lists. Note that the grouping is on the _local_time_ recorded in each file to
    # ensure that these are midnight to midnight means in each location rather than
    # using the UTC times.

    with xarray.open_dataset(each_file) as dat:
        # Reduce to the 2001-2010 focal period along the time axis
        dat = dat.isel(
            time=np.logical_and(
                dat.local_time >= np.datetime64("2001-01-01"),
                dat.local_time <= np.datetime64("2010-12-31"),
            )
        )

        # Daily means using local time for the standard and subdaily gpp and convert
        # from µg C m-2 s-1 to gC m-2 day
        daily_standard_pmod = (
            dat.standard_gpp.groupby("local_time.date").mean() * (60 * 60 * 24) / 1e6
        )
        daily_subdaily_pmod = (
            dat.standard_gpp.groupby("local_time.date").mean() * (60 * 60 * 24) / 1e6
        )

        # Now need to multiply the subdaily value by the soil beta penalty
        daily_subdaily_pmod_with_beta = (
            daily_subdaily_pmod
            * soil_beta.sel(daily_subdaily_pmod.coords).soilmstress_mengoli
        ).compute()

        # Stores for processed data
        standard_daily.append(daily_standard_pmod)
        subdaily_daily.append(daily_subdaily_pmod)
        subdaily_with_beta_daily.append(daily_subdaily_pmod_with_beta)

        print("Processed file: {each_file}")


# Compile and export
standard_daily = xarray.merge(standard_daily)
subdaily_daily = xarray.merge(subdaily_daily)
subdaily_with_beta_daily = xarray.merge(subdaily_with_beta_daily)


dataset = xarray.Dataset(
    standard_daily=standard_daily,
    subdaily_daily=subdaily_daily,
    subdaily_with_beta_daily=subdaily_with_beta_daily,
)

# Export daily GPP grids data by year

out_dir = Path(
    "/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"
    "global_subdaily_gpp/daily_gpp_grids"
)

years, datasets = zip(*dataset.groupby("time.year"))
paths = [out_dir / f"daily_gpp{y}.nc" for y in years]
xarray.save_mfdataset(datasets, paths)
