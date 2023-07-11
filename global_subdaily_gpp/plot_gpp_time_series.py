"""Plot a time series from a global_subdaily_gpp results file
"""
import argparse
from pathlib import Path

import xarray
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np


def plot_gpp_time_series(file: Path, lat_idx: int, lon_idx: int, outfile: Path):
    # Open the provided results file
    data = xarray.open_dataset(file)

    # Get the time series for the specified lat and lon idx
    standard_gpp = data["standard_gpp"][:, lat_idx, lon_idx]
    subdaily_gpp = data["subdaily_gpp"][:, lat_idx, lon_idx]

    # Get the time data
    datetimes = data["time"]

    with PdfPages("multipage_pdf.pdf") as pdf:
        for year in np.unique(datetimes.dt.year):
            fig, axes = plt.subplots(
                ncols=1,
                nrows=12,
                sharex=True,
                sharey=True,
                figsize=(8.27, 11.69),
                layout="constrained",
            )

            for month_idx in np.arange(12):
                # Get indices of data to plot - month within year
                to_plot = (datetimes.dt.year == year) & (
                    datetimes.dt.month == (month_idx + 1)
                )

                this_ax = axes[month_idx]
                this_ax.plot(datetimes.dt.day[to_plot], standard_gpp[to_plot])
                this_ax.plot(datetimes.dt.day[to_plot], subdaily_gpp[to_plot])
                this_ax.text(
                    0.95, 0.95, month_idx + 1, transform=this_ax.transAxes, va="top"
                )

            # Add a titles and axis labels, save figure and close page.
            # fig.tight_layout()
            fig.supylabel("Day of month")
            fig.supxlabel("GPP")
            fig.suptitle(year)
            pdf.savefig()
            plt.close(fig)


if __name__ == "__main__":
    # Parse the arguments provided to main and execute the function.
    parser = argparse.ArgumentParser("plot_gpp_time_series")
    parser.add_argument("--file", type=Path)
    parser.add_argument("--outfile", type=Path)
    parser.add_argument("--lat_idx", type=int)
    parser.add_argument("--lon_idx", type=int)
    args = parser.parse_args()

    plot_gpp_time_series(
        file=args.file,
        outfile=args.outfile,
        lat_idx=args.lat_idx,
        lon_idx=args.lon_idx,
    )
