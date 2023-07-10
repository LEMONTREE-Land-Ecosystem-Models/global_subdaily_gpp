"""This script converts Giulia's details for the 67 Fluxnet sites into an elev input 
file like the ASurf_land_cells.nc required for running global SPLASH
"""
import pandas
import numpy as np

sites = pandas.read_csv("67Sites_geo_info.csv")

# round lat and lon to nearest half degree resolution, which is used as the cell centre
# in the CRU datasets
ll_vals = np.arange(-179.75, 179.75, 0.5)
sites["lat"] = np.array([ll_vals[np.argmin(np.abs(a - ll_vals))] for a in sites.LAT])
sites["lon"] = np.array([ll_vals[np.argmin(np.abs(a - ll_vals))] for a in sites.LON])
sites["cell_id"] = np.arange(1, 68)

sites = sites.drop(columns=["LAT", "LON", "Unnamed: 4"])

# Get the simple xarray indexed by cell_id
sites_indexed = sites.drop(columns=["lat", "lon", "SITE_ID"])
sites_indexed = sites_indexed.set_index("cell_id")
sites_xr = sites_indexed.to_xarray()

# Add the additional coordinates along the cell_id dimension
sites_xr = sites_xr.assign_coords(
    lat=("cell_id", sites.lat),
    lon=("cell_id", sites.lon),
    name=("cell_id", sites.SITE_ID),
)

sites_xr.to_netcdf("67Sites_geo_info.nc")
