library(ncdf4)
library(lubridate)

# Open and read the data
nc_data <- nc_open("67Sites_results.nc")
nc_geo <- nc_open("67Sites_geo_info.nc")

# The geo file is a simple dataframe of one dimensional vectors
nc_geo_df <- data.frame(
    cell_id = ncvar_get(nc_geo, "cell_id"),
    elv = ncvar_get(nc_geo, "elv"),
    lat = ncvar_get(nc_geo, "lat"),
    lon = ncvar_get(nc_geo, "lon"),
    name = ncvar_get(nc_geo, "name")
)

# Get the timeseries data
time_sequence <- seq(
    as.POSIXct("2001-01-01"),
    to = as.POSIXct("2020-12-31"),
    by = "1 day"
)
year <- as.POSIXlt(time_sequence)$year + 1900
pet <- ncvar_get(nc_data, "pet")
aet <- ncvar_get(nc_data, "aet")
pre <- ncvar_get(nc_data, "pre")
cell_id <- ncvar_get(nc_data, "cell_id")

pdf("AI_comparison.pdf", height = 11, width = 8)

par(mfrow = c(5, 2), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))

# Iterate over the sites
for (site_idx in seq_len(nrow(nc_geo_df))) {
    # Site details
    this_site <- as.list(nc_geo_df[site_idx, ])

    # Original result for comparison
    giulia_result <- read.csv(
        file.path(
            "giulia_splash_sites",
            paste0(this_site$name, "_splashv1_2001-2020dfSPLASH_1.csv")
        )
    )

    # Extract the site results from the new data
    new_results <- data.frame(
        aet = aet[, this_site$cell_id],
        pet = pet[, this_site$cell_id],
        pre = pre[, this_site$cell_id]
    )

    # Get annual aridity index for the original and new results
    giulia_ai <- aggregate(
        cbind(total_prec = giulia_result$Prec, total_pet = giulia_result$ep_n),
        by = list(year = year), FUN = sum
    )
    giulia_ai$AI <- giulia_ai$total_pet / giulia_ai$total_prec

    new_ai <- aggregate(
        cbind(total_prec = new_results$pre, total_pet = new_results$pet),
        by = list(year = year), FUN = sum
    )
    new_ai$AI <- new_ai$total_pet / new_ai$total_prec

    # Plot the 20 year AI time series
    plot(AI ~ year,
        data = giulia_ai, type = "s", col = "blue",
        ylim = range(c(giulia_ai$AI, new_ai$AI))
    )
    lines(AI ~ year, data = new_ai, type = "s", col = "red")
    title(main = sprintf(
        "%s: New AI = %0.3f, Giulia AI = %0.3f",
        this_site$name, mean(new_ai$AI), mean(giulia_ai$AI)
    ))
}

dev.off()
