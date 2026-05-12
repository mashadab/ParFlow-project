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
# 4. Download/extract SSEBop actual ET using PyGeoHydro
# --------------------------------------------------
dates = ("2003-10-01", "2004-10-01")

eta = gh.ssebopeta_bygeom(
    bbox_lcc,
    dates=dates,
    geo_crs=conus2_crs
)

print("\nSSEBop ET DataArray:")
print(eta)

# --------------------------------------------------
# 5. Save full gridded daily ET
# --------------------------------------------------
eta.to_netcdf("ssebop_eta_conus2_bbox_daily.nc")

# --------------------------------------------------
# 6. Regional daily mean ET
# --------------------------------------------------
daily_mean = eta.mean(dim=("x", "y"), skipna=True)

df_daily = daily_mean.to_dataframe(name="ET_mm_day").reset_index()
df_daily.to_csv("ssebop_eta_conus2_bbox_daily_mean.csv", index=False)

# --------------------------------------------------
# 7. Monthly total ET
# --------------------------------------------------
monthly_total = daily_mean.resample(time="MS").sum()

df_monthly = monthly_total.to_dataframe(name="ET_mm_month").reset_index()
df_monthly.to_csv("ssebop_eta_conus2_bbox_monthly_total.csv", index=False)

print("\nDaily mean ET:")
print(df_daily.head())

print("\nMonthly total ET:")
print(df_monthly)

print("\nSaved files:")
print("ssebop_eta_conus2_bbox_daily.nc")
print("ssebop_eta_conus2_bbox_daily_mean.csv")
print("ssebop_eta_conus2_bbox_monthly_total.csv")