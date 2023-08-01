# Fitting subdaily global models

This repository contains code to fit the P Model using slow responses on global datasets
using `pyrealm` and to calculate global soil moisture and aridity inputs using SPLASH
v1.0 to calculate the Mengoli et al soil moisture penalty factors for GPP.

## SPLASH calculations of aridity and soil moisture

The Python code version of the SPLASH v1.0 release is used to calculate daily soil
moisture and climatological aridity indices at 0.5° spatial resolution. These values are
calculated using CRU TS4.06 data, provided as 0.5° resolution monthly data:

* Monthly air temperatures are filled to give identical daily values.
* Monthly precipitation is evenly divided across the days of each month.
* Monthly sunshine fraction is calculated as the inverse of CRU cloud fraction
  proportions.

The `compile_cru.py` Python script is used to read the original decadal files for
2000-2019 and reduce the dataset from lat/lon grids to a dataset of landcells identified
by grid ID for processing. This is because the original SPLASH implementation can only
calculate results for a single time series, so iterating over global datasets is much
faster if reduced to only valid inputs.

### Implementation of SPLASH calculation

The SPLASH calculations use a downloaded copy of the original release of the Python
version of SPLASH 1.0. This directory of code files is then called using
`splash_CRU_cli.py` which provides a proper command line interface to call the SPLASH
functionality. The `splash_cru_parallel.pbs.sh` script is then used to submit an array
job to the HPC - each subjob tackles a subset of the global land cells, and the Python
script uses multi-core processing within each subjob to further accelerate the
calculations.

### Validation of SPLASH data

The data and scripts in the `sites_comparison` directory are used to run the SPLASH
model for 67 sites used by Giulia in FluxNet site analyses by running the
`sites_comparison/splash_cru_parallel_67sites.pbs.sh` script on the HPC. This requires
the input file `67Sites_geo_info.csv` containing the details for the 67 sites, created
using `67Sites_to_cell_list.py` and creates the output files `67Sites_results.nc`.

These can then be compared against the outputs of the earlier analyses run by Giulia,
saved in the `giulia_splash_sites` directory. The `AI_comparison.R` file generates site
by site plots showing the annual sequences of aridity index and the climatological 20
year AI for the sites.

### Running global SPLASH calculations

The `splash_cru_parallel.pbs.sh` file is used to run the full set of global 0.5° land
cells, detailed in the `data/Asurf_land_cells.nc` file created by `compile_cru.py`,
through the SPLASH calculations on the HPC.

### Compiling output aridity and soil moisture files

The `compile_aridity_and_soil_moisture.py` and
`compile_aridity_and_soil_moisture.pbs.sh` files are then used to take the resulting
land cells results and calculate global gridded aridity indices and soil moisture time
series.

## P Model fitting

The models are fitted at a 0.5° global resolution to hourly data between 2000-01-01 and
2019-12-31, giving total model dimensions:

* 175320 hourly time points - but see local time notes below.
* 360 latitudinal bands
* 720 longitudinal bands

The models use the following data sources:

* Air temperature from hourly WFDE5 v2 `Tair` data at 0.5° resolution.
* Atmospheric pressure calculated from elevation from WFDE5 v2 `ASurf` at 0.5°
  resolution. These values are constant through time for sites; WFDE5 v2 does provide an
  hourly time series for atmospheric pressure but that is not currently used.
* Vapour pressure deficit calculated from the hourly WFDE5 v2 `Qair` data for specific
  humidity at 0.5° resolution and converted using temperature and pressure as above.
* CO2 concentration data taken from [NOAA monthly global sea level
  data](https://gml.noaa.gov/ccgg/trends/gl_data.html), with data filled to identical
  hourly values within months and across all cells.
* PPFD calculated from the hourly WFDE5 v2 `SWdown` data at 0.5° resolution.
* FAPAR calculated from the daily SNU 0.05° fAPAR data, downsampled using `cdo remapcon`
  to the 0.5° resolution abd then filled to give identical hourly fAPAR values within
  days. Note that some care is need with the downsampling to match the origin and
  orientation of the data to the WFDE5 datasets - the grid of the downsampled data needs
  to be set to match, not just use the default from using `r360x720`.

Gaps in the data inputs are filled by simple forward filling of values along the time
dimension.

The standard and subdaily P Models are fitted to the data:

* Standard P Model (instantaneous acclimation) using `kphio` of 1/8 and the standard C3
  calculations for optimal $\chi$.
* Subdaily P Model (slow acclimation of Jmax, Vcmax and $\xi$) again using `kphio` of
  1/8 and a memory effect for slow responses of $\alpha = 1/15$. The slowly responding
  parameters are set to track noon conditions.

### Local times

The models are divided into longitudinal bands to be run in parallel on the HPC. Because
the input datasets use UTC times, the models are fitted to each column of longitudinal
data independently, using the appropriate local time to correctly identify noon
conditions within the datasets.

Handling this has a couple of effects:

* The local times for the study period do not cover a set of whole days, which at the
  moment is required in the inputs for calculating realised slow responses. The datasets
  are trimmed to 2000-01-02 and 2019-12-30 in order to get the same number of complete
  days using the local calendar.

* The timescale in the output files gives the **local time** for each cell.

### Implementation of P Model calculations

The modelling is implemented using the `global_subdaily_gpp.py` Python script using
`xarray` for data loading and handling and `pyrealm` for running the P Model
calculations. The `global_subdaily_gpp.pbs.sh` shell script is used to run an array job
on the Imperial HPC, allowing different groups of longitudinal bands to be run on
different nodes in parallel. See the documentation in the Python script for details of
the options used in running these jobs.

### Validation

The following commands are used to run the global models for the individual longitudinal
bands containing test data for 7 Fluxnet sites where subdaily models have been run using
Giulia's original code. The odd job array specification is simply being used to submit
an array job with only one subjob to run the specific longitudinal bands containing
these test sites.

``` sh
qsub -J 163-164:2 -v WRITE_PMODEL_INPUTS=TRUE,N_LON_SLICES=720 global_subdaily_models.pbs.sh
qsub -J 190-191:2 -v WRITE_PMODEL_INPUTS=TRUE,N_LON_SLICES=720 global_subdaily_models.pbs.sh
qsub -J 641-642:2 -v WRITE_PMODEL_INPUTS=TRUE,N_LON_SLICES=720 global_subdaily_models.pbs.sh
qsub -J 371-372:2 -v WRITE_PMODEL_INPUTS=TRUE,N_LON_SLICES=720 global_subdaily_models.pbs.sh
qsub -J 181-182:2 -v WRITE_PMODEL_INPUTS=TRUE,N_LON_SLICES=720 global_subdaily_models.pbs.sh
qsub -J 187-188:2 -v WRITE_PMODEL_INPUTS=TRUE,N_LON_SLICES=720 global_subdaily_models.pbs.sh
qsub -J 408-409:2 -v WRITE_PMODEL_INPUTS=TRUE,N_LON_SLICES=720 global_subdaily_models.pbs.sh
```

### Running global calculations

The complete set of longitudinal bands can then be run using:

``` sh
qsub -J 0-29 -v N_LON_SLICES=30 global_subdaily_models.pbs.sh
```

This runs 30 nodes, each tackling a set of 24 longitudinal bands. Each band is written
out as a separate NetCDF file.

### Compiling output files

The `daily_gpp_grids` directory contains compiled GPP estimates in gC m-2 day from
2001 - 2010. The values include the predictions of the standard P Model GPP, the
subdaily P Model GPP and the subdaily P Model GPP applying the Mengoli et al soil beta
correction.

The hourly GPP output data are processed as follows:

* The hourly outputs are extracted from each longitudinal band for the focal time
  period.
* The hourly values are grouped into daily sets of observations using the `local_time`
  coordinates in those files. This ensures that all grid cells are true midnight to
  midnight means, rather than UTC midnight to midnight, regardless of their longitude.
* The mean GPP within days is calculated for both the standard and subdaily P Model
  outputs and converted from µgC m-2 s-1 to gC m-2 day-1
  ($(60 \times 60 \times 24) / 1e6$).
* The matching Mengoli soil beta factor is the extracted for the longitudinal band and
  time period and applied to the subdaily estimate of daily GPP.

The resulting data from longitudinal bands are then compiled into global grids and saved
into annual files containing each of the three GPP estimates.

For some reason, the `daily_gpp_compiler.py` script has a huge RAM requirement. The
eventual annual files are only 2.2 GB, but each annual subjob requires ~540 GB of RAM.
This seems _very_ odd to me.
