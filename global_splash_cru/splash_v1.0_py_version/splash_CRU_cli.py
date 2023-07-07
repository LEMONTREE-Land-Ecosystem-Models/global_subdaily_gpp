"""Run SPLASH v1 on global grid data

# Single cell run for testing (cell is ~ Crediton)
python splash_CRU_cli.py \
    --tmp ../../data/cru_2001_2020_tmp.nc \
    --pre ../../data/cru_2001_2020_pre.nc \
    --cld ../../data/cru_2001_2020_cld.nc \
    --elv ../../data/Asurf_land_cells.nc \
    --cell_ids 202672

# Select a subset block of cells to run
python splash_CRU_cli.py \
    --tmp ../../data/cru_2001_2020_tmp.nc \
    --pre ../../data/cru_2001_2020_pre.nc \
    --cld ../../data/cru_2001_2020_cld.nc \
    --elv ../../data/Asurf_land_cells.nc \
    --block_number 0
    --total_blocks 30

"""

import argparse
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
from itertools import repeat


from multiprocess.pool import Pool
import xarray
import numpy as np

from splash import SPLASH


@dataclass
class SpinUpData:
    """Lightweight drop in replacement for the original SPLASH DATA class, which has
    hard coded paths and all sorts of other difficulties. This class is just used to
    package up one year of data to pass to the SPLASH.spinup method."""

    pn_vec: np.ndarray
    tair_vec: np.ndarray
    sf_vec: np.ndarray
    year: int
    npoints: int = field(init=False)
    num_lines: int = field(init=False)

    def __post_init__(self):
        self.npoints = len(self.pn_vec)
        self.num_lines = len(self.pn_vec)


def run_one_cell(cell, climate_data, dates):
    """Run SPLASH on a single cell.

    This function takes a cell object (with lat, lon and cell_id information) and then
    uses that information to extract the monthly data from the input climate data, pad
    it to daily and then run the spinup and calculate the time series.
    """
    cell_data = {}

    # Extract the monthly data and fill to daily
    for ky, ds in climate_data.items():
        this_var = ds.sel(lat=float(cell.lat), lon=float(cell.lon))
        cell_data[ky] = this_var.resample(time="D").ffill()

    # Calculate input variables and reduce to numpy.
    # Divide precipitation by days in month
    pre = cell_data["pre"].values / dates.dt.days_in_month.values
    # Calculate sunshine fraction
    sf = 1 - (cell_data["cld"].values / 100)
    tmp = cell_data["tmp"].values

    # Create the drop in data replacement for the first year
    in_year_one = dates < np.datetime64("2002-01-01")
    spin_up_data = SpinUpData(
        tair_vec=tmp[in_year_one],
        pn_vec=pre[in_year_one],
        sf_vec=sf[in_year_one],
        year=2001,
    )

    # Initialise SPLASH object
    splash = SPLASH(lat=float(cell.lat), elv=float(cell.values))

    # Run the spin up year
    splash.spin_up(spin_up_data)
    curr_wn = splash.wn

    # Loop through the full time sequence
    dayofyear = dates.dt.dayofyear
    year = dates.dt.year

    # Create stores for the output
    cond_out = np.full_like(pre, np.nan)
    eet_out = np.full_like(pre, np.nan)
    pet_out = np.full_like(pre, np.nan)
    aet_out = np.full_like(pre, np.nan)
    wn_out = np.full_like(pre, np.nan)
    ro_out = np.full_like(pre, np.nan)

    for idx, (doy, yr) in enumerate(zip(dayofyear, year)):
        # Calculate values and store
        splash.run_one_day(
            n=int(doy),
            y=int(yr),
            wn=curr_wn,
            sf=sf[idx],
            tc=tmp[idx],
            pn=pre[idx],
        )

        cond_out[idx] = splash.cond
        eet_out[idx] = splash.eet
        pet_out[idx] = splash.pet
        aet_out[idx] = splash.aet
        wn_out[idx] = splash.wn
        ro_out[idx] = splash.ro

        # update the wn for the following day
        curr_wn = splash.wn

    print(f"Cell: {cell.cell_id.values} completed")

    # Return values of interest- also include input variables to make it easier to
    # calculate aridity index etc.
    return dict(
        cell=cell,
        cond=cond_out,
        eet=eet_out,
        pet=pet_out,
        aet=aet_out,
        wn=wn_out,
        ro=ro_out,
        pre=pre,
        sf=sf,
        tmp=tmp,
    )


def splash_cru_cli(
    tmp: Path,
    pre: Path,
    cld: Path,
    elv: Path,
    output_file: Path,
    block_number: Optional[int] = 1,
    total_blocks: Optional[int] = 30,
    cell_ids: Optional[tuple[int]] = None,
    n_process: Optional[int] = None,
):
    """Implement SPLASH calculations from monthly CRU data.

    The function requires paths to the tmp, pre and cld CRU grids and to an elevation
    file containing the elevations of land cells, indexed by cell_id, lat and lon. A
    specified block of the land cells can be run, or a tuple of specific cell ids can
    also be run. The n_process argument can be supplied to apply multiprocessing of the
    resulting set of cell_ids.

    """

    # Load elevation and latitude data for land cells
    elv_data = xarray.load_dataarray(elv)

    if cell_ids is None:
        # Create block id numbers starting at 1:total_blocks for compatibility with PBS
        # job array values, truncated to the actual total number of cells
        n_cells = len(elv_data)
        cells_per_block = int(np.ceil(n_cells / total_blocks))
        block_ids = np.repeat(np.arange(1, total_blocks + 1), cells_per_block)
        block_ids = block_ids[0:n_cells]

        # Subset to the cells being handled in this block
        elv_data = elv_data.loc[block_ids == block_number]
    else:
        # Alternatively run a specific set of cells
        elv_data = elv_data.sel(cell_id=list(cell_ids))

    # Load the climate source files
    climate_sources = {"pre": pre, "tmp": tmp, "cld": cld}
    climate_data = {}

    for ky, src in climate_sources.items():
        raw_data = xarray.load_dataarray(src)

        # reset the time indices from mid month coords in data files to start month
        # values (and then back to nano-second precision to suppress a warning)
        month_start = raw_data.time.values.astype("datetime64[M]").astype(
            "datetime64[ns]"
        )
        raw_data.coords["time"] = month_start

        # Copy the start of the last month to the end last month to make forward fill to
        # last day
        pad_data = raw_data.isel(time=-1)
        pad_data.coords["time"] = np.datetime64("2020-12-31").astype("datetime64[ns]")
        climate_data[ky] = xarray.concat([raw_data, pad_data], dim="time")

    # Get the dates as a DataArray to use the datetime accessor functions
    dates = xarray.DataArray(
        np.arange(month_start[0], np.datetime64("2021-01-01"), np.timedelta64("1", "D"))
        .astype("datetime64[D]")
        .astype("datetime64[ns]")
    )

    # Iterate over the subset of cells
    if n_process is not None:
        # Use a pool of processes to run the set of cells
        with Pool(n_process) as pool:
            result = pool.starmap_async(
                run_one_cell, zip(elv_data, repeat(climate_data), repeat(dates))
            )

            results = result.get()

    else:
        # Run the cells in series
        results = []
        for cell in elv_data:
            cell_results = run_one_cell(
                cell=cell, climate_data=climate_data, dates=dates
            )
            results.append(cell_results)

    # Build results into a output file
    cell_id = [int(r["cell"].cell_id) for r in results]
    # lat = [float(r["cell"].lat) for r in results]
    # lon = [float(r["cell"].lon) for r in results]

    out_data = {}

    # Compile an xarray of each of the main output variables, indexed by cell_id and
    # time - also include precipitation to make it easier to calculate aridity index
    for var in ("cond", "eet", "pet", "aet", "wn", "ro", "pre", "sf", "tmp"):
        this_data_var = [r[var] for r in results]
        this_data_var = np.stack(this_data_var)
        out_data[var] = xarray.DataArray(
            this_data_var.astype(np.float32),
            coords={"cell_id": cell_id, "time": dates.values},
        )

    out_data_ds = xarray.Dataset(out_data)
    out_data_ds.to_netcdf(output_file)

    print("Cells completed")


if __name__ == "__main__":
    # Parse the arguments provided to main and execute the function.
    parser = argparse.ArgumentParser("splash_cli")
    parser.add_argument("--tmp", type=Path)
    parser.add_argument("--cld", type=Path)
    parser.add_argument("--pre", type=Path)
    parser.add_argument("--elv", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--block_number", type=int, default=None)
    parser.add_argument("--total_blocks", type=int, default=None)
    parser.add_argument("--cell_ids", type=int, nargs="+", default=None)
    parser.add_argument("--n_process", type=int, default=None)
    args = parser.parse_args()

    splash_cru_cli(
        tmp=args.tmp,
        pre=args.pre,
        cld=args.cld,
        elv=args.elv,
        output_file=args.out,
        block_number=args.block_number,
        total_blocks=args.total_blocks,
        cell_ids=args.cell_ids,
        n_process=args.n_process,
    )
