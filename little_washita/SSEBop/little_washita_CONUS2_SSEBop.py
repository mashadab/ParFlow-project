# extract_ssebop_conus2_bbox.py

import pygeohydro as gh
import pandas as pd
from pyproj import CRS, Transformer
from shapely.geometry import box

# --------------------------------------------------
# 1. Input lat/lon bounding box
# Format: [[south_lat, west_lon], [north_lat, east_lon]]
# --------------------------------------------------
latlon_bounds = [[34.739932, -98.420500], [35.031552, -97.721777]]

# --------------------------------------------------
# 2. CONUS2 Lambert Conformal Conic CRS
# --------------------------------------------------
conus2_crs = CRS.from_proj4(
    "+proj=lcc "
    "+lat_1=30 +lat_2=60 "
    "+lon_0=-97.0 "
    "+lat_0=40.0000076294444 "
    "+a=6370000.0 +b=6370000 "
    "+units=m +no_defs"
)

wgs84 = CRS.from_epsg(4326)

'''
xmin = -127500.343
ymin = -570338.692
xmax = -64500.312
ymax = -537664.372
'''


# --------------------------------------------------
# 3. Convert lat/lon bbox to CONUS2 LCC bbox
# --------------------------------------------------
transformer = Transformer.from_crs(wgs84, conus2_crs, always_xy=True)

south_lat, west_lon = latlon_bounds[0]
north_lat, east_lon = latlon_bounds[1]

corners_lonlat = [
    (west_lon, south_lat),  # southwest
    (west_lon, north_lat),  # northwest
    (east_lon, north_lat),  # northeast
    (east_lon, south_lat),  # southeast
]

corners_lcc = [transformer.transform(lon, lat) for lon, lat in corners_lonlat]

xs = [x for x, y in corners_lcc]
ys = [y for x, y in corners_lcc]

xmin, xmax = min(xs), max(xs)
ymin, ymax = min(ys), max(ys)

bbox_lcc = box(xmin, ymin, xmax, ymax)

print("CONUS2 LCC corners:")
for i, pt in enumerate(corners_lcc, start=1):
    print(f"corner {i}: x={pt[0]:.3f}, y={pt[1]:.3f}")

print("\nCONUS2 LCC bounding box:")
print(f"xmin = {xmin:.3f}")
print(f"ymin = {ymin:.3f}")
print(f"xmax = {xmax:.3f}")
print(f"ymax = {ymax:.3f}")

# --------------------------------------------------
# Bounding box dimensions
# --------------------------------------------------
width_m = xmax - xmin
height_m = ymax - ymin

width_km = width_m / 1000
height_km = height_m / 1000

print("\nBounding box dimensions:")
print(f"Width  = {width_m:.2f} m ({width_km:.2f} km)")
print(f"Height = {height_m:.2f} m ({height_km:.2f} km)")
print(f"Area   = {(width_m * height_m) / 1e6:.2f} km²")

# --------------------------------------------------
# 4. Download/extract SSEBop actual ET using PyGeoHydro
# --------------------------------------------------
'''dates = ("2003-10-01", "2004-10-01")

eta = gh.ssebopeta_bygeom(
    bbox_lcc,
    dates=dates,
    geo_crs=conus2_crs
)

print("\nSSEBop ET DataArray:")
print(eta)
'''
# --------------------------------------------------
# 4. Download/extract SSEBop actual ET
#    Daily reading with bad TIFF skipping
# --------------------------------------------------
import xarray as xr
import numpy as np

start_date = "2003-10-01"
end_date   = "2004-10-01"

days = pd.date_range(
    start=start_date,
    end=end_date,
    freq="D",
    inclusive="left"
)

arrays = []
failed_days = []

for day in days:

    d0 = day.strftime("%Y-%m-%d")
    d1 = (day + pd.Timedelta(days=1)).strftime("%Y-%m-%d")

    try:
        print(f"Reading {d0}")

        eta_day = gh.ssebopeta_bygeom(
            bbox_lcc,
            dates=(d0, d1),
            geo_crs=conus2_crs
        )

        # Skip empty or invalid datasets
        if eta_day is None:
            print(f"Skipping {d0}: returned None")
            failed_days.append(d0)
            continue

        # Make sure it actually has values
        if eta_day.size == 0:
            print(f"Skipping {d0}: empty array")
            failed_days.append(d0)
            continue

        # Skip all-NaN rasters
        vals = eta_day.values
        if np.isnan(vals).all():
            print(f"Skipping {d0}: all NaNs")
            failed_days.append(d0)
            continue

        arrays.append(eta_day)

    except Exception as e:
        # catches RasterioIOError and bad TIFFs
        print(f"Skipping {d0}: {str(e)}")
        failed_days.append(d0)
        continue

# --------------------------------------------------
# Combine valid days
# --------------------------------------------------
if len(arrays) == 0:
    raise RuntimeError(
        "No valid SSEBop rasters were downloaded."
    )

eta = xr.concat(arrays, dim="time")

print("\nFinal dataset:")
print(eta)

print(f"\nValid days: {len(arrays)}")
print(f"Skipped days: {len(failed_days)}")

# --------------------------------------------------
# 5. Save gridded daily ET
# --------------------------------------------------
eta.to_netcdf("ssebop_eta_conus2_bbox_daily.nc")

# --------------------------------------------------
# 6. Regional daily mean ET
# --------------------------------------------------
daily_mean = eta.mean(
    dim=("x", "y"),
    skipna=True
)

df_daily = (
    daily_mean
    .to_dataframe(name="ET_mm_day")
    .reset_index()
)

df_daily.to_csv(
    "ssebop_eta_conus2_bbox_daily_mean.csv",
    index=False
)

# --------------------------------------------------
# 66666666. Regional daily sum ET  #rerun this block
# --------------------------------------------------
daily_sum = eta.sum(
    dim=("x", "y"),
    skipna=True
) * 1e3 * 1e3  #km to m conversion  m3/day

df_daily_sum = (
    daily_sum
    .to_dataframe(name="ET_m3_day")
    .reset_index()
)

df_daily_sum.to_csv(
    "ssebop_eta_conus2_bbox_daily_sum.csv",
    index=False
)


# Keep only unique dates
df_daily_sum_unique = (
    df_daily_sum
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

df_daily_sum_unique.to_csv(
    "ssebop_eta_conus2_bbox_daily_sum_unique.csv",
    index=False
)

print(f"Unique dates: {len(df_daily_sum)}")


# Keep only unique dates
daily_sum_mpday = eta.sum(
    dim=("x", "y"),
    skipna=True
)/(width_km*height_km)  #km to m conversion  m/day

df_daily_sum_mpday = (
    daily_sum_mpday
    .to_dataframe(name="ET_m_day")
    .reset_index()
)

# Keep only unique dates
df_daily_sum_unique_mpday = (
    df_daily_sum_mpday
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)


df_daily_sum_unique_mpday.to_csv(
    "ssebop_eta_conus2_bbox_daily_sum_unique_mpday.csv",
    index=False
)


# --------------------------------------------------
# 7. Monthly total ET
# --------------------------------------------------
monthly_total = (
    daily_mean
    .resample(time="MS")
    .sum()
)

df_monthly = (
    monthly_total
    .to_dataframe(name="ET_mm_month")
    .reset_index()
)

df_monthly.to_csv(
    "ssebop_eta_conus2_bbox_monthly_total.csv",
    index=False
)

# --------------------------------------------------
# 8. Save skipped days log
# --------------------------------------------------
pd.DataFrame(
    {"failed_days": failed_days}
).to_csv(
    "ssebop_failed_days.csv",
    index=False
)

print("\nSaved files:")
print("ssebop_eta_conus2_bbox_daily.nc")
print("ssebop_eta_conus2_bbox_daily_mean.csv")
print("ssebop_eta_conus2_bbox_monthly_total.csv")
print("ssebop_failed_days.csv")