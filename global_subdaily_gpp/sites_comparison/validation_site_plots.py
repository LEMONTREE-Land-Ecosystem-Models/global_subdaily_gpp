from pathlib import Path

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


pdf = PdfPages("temp.pdf")

for sitename, fpath in giulia_sites.items():
    # Load the data and get datetimes and decimal days in month
    local_data = pandas.read_csv(fpath)
    local_datetime = pandas.to_datetime(local_data["TIME"])
    local_decimal_day_in_month = (
        local_datetime.dt.day
        + local_datetime.dt.hour / 24
        + local_datetime.dt.minute / 24 * 60
    )

    # Get site info
    site_info = sites.loc[sitename]

    print(site_info)
    print(
        local_datetime.min(),
        local_datetime.max(),
        np.diff(local_datetime[0:2]).astype("timedelta64[m]"),
    )

    # Open the provided results file and subset to the time series and cell required
    global_data = xarray.open_dataset(
        f"../results/gpp_data_{int(site_info.lon_idx.values)}.nc"
    )

    global_subset = global_data.sel(
        time=slice(local_datetime.min().to_numpy(), local_datetime.max().to_numpy()),
        lat=site_info.lat_band,
        lon=site_info.lon_band,
    )

    # Get the global time series and decimal day in month to plot by month
    global_datetime = global_subset["time"]
    global_decimal_day_in_month = (
        global_datetime.dt.day
        + global_datetime.dt.hour / 24
        + global_datetime.dt.minute / 24 * 60
    )

    for year in np.unique(local_datetime.dt.year):
        fig, axes = plt.subplots(
            ncols=1,
            nrows=12,
            sharex=True,
            sharey=True,
            figsize=(8.27, 11.69),
            layout="constrained",
        )

        for month_idx in np.arange(12):
            # Get indices of data to plot - month within year - for both local and
            # global time series.
            local_to_plot = (local_datetime.dt.year == year) & (
                local_datetime.dt.month == (month_idx + 1)
            )
            global_to_plot = (global_datetime.dt.year == year) & (
                global_datetime.dt.month == (month_idx + 1)
            )

            # Plot the predictions
            this_ax = axes[month_idx]
            this_ax.plot(
                local_decimal_day_in_month[local_to_plot],
                local_data["GPPp"][local_to_plot],
                label="Standard P Model - site",
            )
            this_ax.plot(
                local_decimal_day_in_month[local_to_plot],
                local_data["GPPpOpt"][local_to_plot],
                label="Standard P Model - site",
            )
            this_ax.plot(
                global_decimal_day_in_month[global_to_plot],
                global_subset["standard_gpp"][global_to_plot],
                label="Standard P Model - site",
            )
            this_ax.plot(
                global_decimal_day_in_month[global_to_plot],
                global_subset["subdaily_gpp"][global_to_plot],
                label="Standard P Model - site",
            )

            # Add a month number and legend
            this_ax.text(
                0.95, 0.95, month_idx + 1, transform=this_ax.transAxes, va="top"
            )
            if month_idx == 0:
                this_ax.legend()

        # Add a titles and axis labels, save figure and close page.
        # fig.tight_layout()
        fig.supxlabel("Decimal day of month")
        fig.supylabel("GPP")
        fig.suptitle(f"{sitename} : {year}")
        pdf.savefig()
        plt.close(fig)

pdf.close()
