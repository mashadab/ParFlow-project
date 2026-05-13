#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 11:06:57 2026

@author: ms6985
"""

# extract_ssebop_latlon_bbox.py

import pygeohydro as gh
import pandas as pd
import numpy as np
import xarray as xr
from pyproj import CRS, Geod
from shapely.geometry import box

# --------------------------------------------------
# 1. Input lat/lon bounding box
# Format: [[south_lat, west_lon], [north_lat, east_lon]]
# --------------------------------------------------
latlon_bounds = [[42.3433620411624929,-72.2925311203326686], [42.692484048100727,-71.9502074688803788]]


south_lat, west_lon = latlon_bounds[0]
north_lat, east_lon = latlon_bounds[1]

# --------------------------------------------------
# 2. Build regular lat/lon bbox
# shapely uses: box(xmin, ymin, xmax, ymax)
# where x = lon, y = lat
# --------------------------------------------------
bbox_latlon = box(west_lon, south_lat, east_lon, north_lat)

wgs84 = CRS.from_epsg(4326)

print("Lat/lon bounding box:")
print(bbox_latlon.bounds)

# --------------------------------------------------
# 3. Estimate bbox dimensions in meters
# --------------------------------------------------
geod = Geod(ellps="WGS84")

mid_lat = 0.5 * (south_lat + north_lat)
mid_lon = 0.5 * (west_lon + east_lon)

_, _, width_m = geod.inv(west_lon, mid_lat, east_lon, mid_lat)
_, _, height_m = geod.inv(mid_lon, south_lat, mid_lon, north_lat)

width_km = width_m / 1000.0
height_km = height_m / 1000.0
area_m2 = width_m * height_m

print("\nApproximate bbox dimensions:")
print(f"Width  = {width_m:.2f} m ({width_km:.2f} km)")
print(f"Height = {height_m:.2f} m ({height_km:.2f} km)")
print(f"Area   = {area_m2 / 1e6:.2f} km²")

# --------------------------------------------------
# 4. Date range
# --------------------------------------------------
start_date = "2019-01-01"
end_date   = "2023-01-01"

days = pd.date_range(
    start=start_date,
    end=end_date,
    freq="D",
    inclusive="left"
)

# --------------------------------------------------
# 5. Download/extract SSEBop daily ET
#    Skip bad/non-TIFF days
# --------------------------------------------------
arrays = []
failed_days = []

for day in days:

    d0 = day.strftime("%Y-%m-%d")
    d1 = (day + pd.Timedelta(days=1)).strftime("%Y-%m-%d")

    try:
        print(f"Reading {d0}")

        eta_day = gh.ssebopeta_bygeom(
            bbox_latlon,
            dates=(d0, d1),
            geo_crs=wgs84
        )

        if eta_day is None:
            print(f"Skipping {d0}: returned None")
            failed_days.append(d0)
            continue

        if eta_day.size == 0:
            print(f"Skipping {d0}: empty array")
            failed_days.append(d0)
            continue

        if np.isnan(eta_day.values).all():
            print(f"Skipping {d0}: all NaNs")
            failed_days.append(d0)
            continue

        arrays.append(eta_day)

    except Exception as e:
        print(f"Skipping {d0}: {str(e)}")
        failed_days.append(d0)
        continue

# --------------------------------------------------
# 6. Combine valid daily rasters
# --------------------------------------------------
if len(arrays) == 0:
    raise RuntimeError("No valid SSEBop rasters were downloaded.")

eta = xr.concat(arrays, dim="time")

# remove repeated dates if any
_, unique_idx = np.unique(
    eta["time"].values,
    return_index=True
)

eta = eta.isel(time=np.sort(unique_idx))

print("\nFinal SSEBop dataset:")
print(eta)

print(f"\nValid days: {eta.sizes['time']}")
print(f"Skipped days: {len(failed_days)}")

# --------------------------------------------------
# 7. Save gridded daily ET
# eta units are mm/day
# --------------------------------------------------
eta.to_netcdf("ssebop_eta_latlon_bbox_daily.nc")

# --------------------------------------------------
# 8. Regional daily mean ET
# mm/day
# --------------------------------------------------
daily_mean_mm = eta.mean(
    dim=("x", "y"),
    skipna=True
)

df_daily_mean = (
    daily_mean_mm
    .to_dataframe(name="ET_mm_day")
    .reset_index()
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

df_daily_mean.to_csv(
    "ssebop_eta_latlon_bbox_daily_mean_mmday.csv",
    index=False
)

# --------------------------------------------------
# 9. Regional daily mean ET
# m/day
# --------------------------------------------------
daily_mean_m = daily_mean_mm / 1000.0

df_daily_mean_m = (
    daily_mean_m
    .to_dataframe(name="ET_m_day")
    .reset_index()
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

df_daily_mean_m.to_csv(
    "ssebop_eta_latlon_bbox_daily_mean_mday.csv",
    index=False
)

# --------------------------------------------------
# 10. Regional daily total ET volume
# m3/day
# --------------------------------------------------
# SSEBop eta is mm/day.
# Approximate cell area is 1 km2 = 1,000,000 m2.
# Volume per cell = ET_mm_day / 1000 * 1,000,000 = ET_mm_day * 1000 m3/day.
cell_area_m2 = 1000.0 * 1000.0

daily_sum_m3 = eta.sum(
    dim=("x", "y"),
    skipna=True
) / 1000.0 * cell_area_m2

df_daily_sum_m3 = (
    daily_sum_m3
    .to_dataframe(name="ET_m3_day")
    .reset_index()
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

df_daily_sum_m3.to_csv(
    "ssebop_eta_latlon_bbox_daily_sum_m3day.csv",
    index=False
)

# --------------------------------------------------
# 11. Monthly total ET
# mm/month and m/month
# --------------------------------------------------
monthly_total_mm = daily_mean_mm.resample(time="MS").sum()
monthly_total_m = monthly_total_mm / 1000.0

df_monthly_mm = (
    monthly_total_mm
    .to_dataframe(name="ET_mm_month")
    .reset_index()
)

df_monthly_m = (
    monthly_total_m
    .to_dataframe(name="ET_m_month")
    .reset_index()
)

df_monthly_mm.to_csv(
    "ssebop_eta_latlon_bbox_monthly_total_mm.csv",
    index=False
)

df_monthly_m.to_csv(
    "ssebop_eta_latlon_bbox_monthly_total_m.csv",
    index=False
)

# --------------------------------------------------
# 12. Save failed days log
# --------------------------------------------------
pd.DataFrame(
    {"failed_days": failed_days}
).to_csv(
    "ssebop_latlon_failed_days.csv",
    index=False
)

print("\nSaved files:")
print("ssebop_eta_latlon_bbox_daily.nc")
print("ssebop_eta_latlon_bbox_daily_mean_mmday.csv")
print("ssebop_eta_latlon_bbox_daily_mean_mday.csv")
print("ssebop_eta_latlon_bbox_daily_sum_m3day.csv")
print("ssebop_eta_latlon_bbox_monthly_total_mm.csv")
print("ssebop_eta_latlon_bbox_monthly_total_m.csv")
print("ssebop_latlon_failed_days.csv")