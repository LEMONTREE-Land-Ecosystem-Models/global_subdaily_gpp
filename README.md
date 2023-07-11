# Fitting subdaily global models

This repository contains code to fit the P Model using slow responses on global datasets
using `pyrealm`.

## P Model fitting

The models are fitted at a 0.5° global resolution to hourly data between 2000-01-01 and
2019-12-31, giving total model dimensions:

* 175320 hourly time points
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
  days.

The standard and subdaily P Models are fitted to the data:

* Standard P Model (instantaneous acclimation) using `kphio` of 1/8 and the standard C3
  calculations for optimal $\chi$.
* Subdaily P Model (slow acclimation of Jmax, Vcmax and $\xi$) again using `kphio` of
  1/8 and a memory effect for slow responses of $\alpha = 1/15$. The slowly responding
  parameters are set to track noon conditions.

The models are divided into longitudinal bands to be run in parallel on the HPC. Because
the input datasets use UTC times, the models are fitted to each column of longitudinal
data independently, using the appropriate local time to correctly identify noon
conditions within the datasets.

## SPLASH calculations of aridity and soil moisture

The Python code version of the SPLASH v1.0 release is used to calculate daily soil
moisture and climatological aridity indices at 0.5° spatial resolution. These values are
calculated using CRU TS4.06 data, provided as 0.5° resolution monthly data:

* Monthly air temperatures are filled to give identical daily values.
* Monthly precipitation is evenly divided across the days of each month.
* Monthly sunshine fraction is calculated as the inverse of CRU cloud fraction
  proportions.

The datasets are reduced to land cells only and then these are processed in parallel
blocks, using multiple cores with each block to further accelerate the calculations.

### Validation of SPLASH data

The data and scripts in the `sites_comparison` directory are used to run the SPLASH
model for 67 sites used by Giulia in FluxNet site analyses. These can then be compared
against the outputs of the earlier analyses.
