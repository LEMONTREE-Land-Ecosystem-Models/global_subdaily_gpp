from pathlib import Path
from itertools import product

import xarray
from matplotlib.backends.backend_pdf import PdfPages
import pandas
import numpy as np
import matplotlib.pyplot as plt

# Get files and extract site name
giulia_sites = list(Path("giulia_data").glob("*.csv"))
giulia_sites = {f.as_posix()[19:25]: f for f in giulia_sites}

# Match to site data
sites = pandas.read_csv("fluxnet_site_geo_info.csv")
sites = sites.iloc[:, 0:4]
sites = sites.set_index(sites.SITE_ID)

# Match lon and lat to nearest 0.5° bands and get the index of the cells along lon and
# lat dimensions
lon_bands = np.arange(-179.75, 180, 0.5)
sites["lon_band"] = np.array(
    [lon_bands[np.argmin(np.abs(a - lon_bands))] for a in sites.LON]
)
sites["lon_idx"] = [int(np.where(lon_bands == a)[0][0]) for a in sites.lon_band]

lat_bands = np.arange(-89.75, 90, 0.5)
sites["lat_band"] = np.array(
    [lat_bands[np.argmin(np.abs(a - lat_bands))] for a in sites.LAT]
)
sites["lat_idx"] = [int(np.where(lat_bands == a)[0][0]) for a in sites.lat_band]


pdf = PdfPages("site_local_global_gpp_comparison.pdf")

for sitename, fpath in giulia_sites.items():
    # Load the data and get datetimes and decimal days in month
    local_data = pandas.read_csv(fpath)
    local_datetime = pandas.to_datetime(local_data["TIME"])
    local_decimal_day_in_month = (
        local_datetime.dt.day
        + local_datetime.dt.hour / 24
        + local_datetime.dt.minute / (24 * 60)
    )

    # Get site info
    site_info = sites.loc[sitename]

    print(site_info)
    print(
        local_datetime.min(),
        local_datetime.max(),
        np.diff(local_datetime[0:2]).astype("timedelta64[m]"),
    )

    # Open the provided results file if found
    global_path = Path(f"../results/gpp_data_{site_info.lon_band}.nc")

    if not global_path.exists():
        print("Missing global results. Skipping.")
        continue

    global_data = xarray.open_dataset(global_path)

    # Subset to the time series and cell required
    global_subset = global_data.sel(
        lat=site_info.lat_band,
        lon=site_info.lon_band,
    )

    global_subset = global_subset.where(
        np.logical_and(
            global_data.local_time > local_datetime.min().to_numpy(),
            global_data.local_time < local_datetime.max().to_numpy(),
        ),
        drop=True,
    )

    # Get the local datetime data from the global time series and decimal day in month
    # to plot by month
    global_datetime = global_subset["local_time"]
    global_decimal_day_in_month = (
        global_datetime.dt.day
        + global_datetime.dt.hour / 24
        + global_datetime.dt.minute / (24 * 60)
    )

    # Plot in pairs of local over global to make it easier to compare data.
    month_blocks = np.split(np.arange(1, 13), 4)
    all_years = np.unique(local_datetime.dt.year)
    month_names = [
        d.item().strftime("%b")
        for d in np.arange(
            np.datetime64("2000-01"), np.datetime64("2001-01"), np.timedelta64(1, "M")
        )
    ]

    for year, months in product(all_years, month_blocks):
        fig, axes = plt.subplots(
            ncols=1,
            nrows=2 * len(months),
            sharex=True,
            sharey=True,
            figsize=(8.27, 11.69),
            layout="constrained",
        )

        for plot_idx, this_month in enumerate(months):
            # Plot the predictions
            print(plot_idx, this_month, 2 * plot_idx, 2 * plot_idx + 1)
            local_ax = axes[2 * plot_idx]
            global_ax = axes[2 * plot_idx + 1]

            # Get indices of data to plot from local data - month within year
            local_to_plot = (local_datetime.dt.year == year) & (
                local_datetime.dt.month == (this_month)
            ).values

            # Plot site based local predictions on local axis
            local_ax.plot(
                local_decimal_day_in_month[local_to_plot],
                local_data["GPPp"][local_to_plot] * 12.011,
                label="Non-acclimating P Model",
                color="grey",
                linewidth=1,
            )
            local_ax.plot(
                local_decimal_day_in_month[local_to_plot],
                local_data["GPPpOpt"][local_to_plot] * 12.011,
                label="Acclimating P Model",
                color="red",
                linewidth=1,
            )

            # Add a plot label
            local_ax.text(
                0.95,
                0.95,
                f"{month_names[this_month - 1]}: Giulia local data",
                transform=local_ax.transAxes,
                va="top",
                ha="right",
            )

            if plot_idx == 0:
                local_ax.legend()

            # Get indices of data to plot from global data
            global_to_plot = (global_datetime.dt.year == year) & (
                global_datetime.dt.month == (this_month)
            ).data

            # Plot global predictions on global axis.
            global_ax.plot(
                global_decimal_day_in_month[global_to_plot],
                global_subset["standard_gpp"][global_to_plot],
                label="Non-acclimating P Model",
                color="grey",
                linewidth=1,
            )
            global_ax.plot(
                global_decimal_day_in_month[global_to_plot],
                global_subset["subdaily_gpp"][global_to_plot],
                label="Acclimating P Model",
                color="red",
                linewidth=1,
            )

            # Add a plot label
            global_ax.text(
                0.95,
                0.95,
                f"{month_names[this_month - 1]}: Pyrealm global data",
                transform=global_ax.transAxes,
                va="top",
                ha="right",
            )

            if plot_idx == 0:
                local_ax.legend()

        # Add a titles and axis labels, save figure and close page.
        # fig.tight_layout()
        fig.supxlabel("Decimal day of month")
        fig.supylabel("GPP (µg C m-2 s-1)")
        fig.suptitle(f"{sitename} : {year}")
        pdf.savefig()
        plt.close(fig)

pdf.close()
