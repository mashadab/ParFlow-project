#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# --------------------------------------------------
# Plot settings
# --------------------------------------------------
plt.rcParams.update({"font.size": 20})
plt.rcParams.update({"font.family": "Serif"})

red   = [190/255,  30/255,  45/255]
blue  = [ 30/255, 144/255, 255/255]
green = [  0/255, 166/255,  81/255]
gray  = [100/255, 100/255, 100/255]

# --------------------------------------------------
# 1. Model domain area
# --------------------------------------------------
model_width_m  = 64000.0
model_height_m = 32000.0
model_area_m2  = model_width_m * model_height_m

# --------------------------------------------------
# 2. Time axis
# --------------------------------------------------
start_datetime = datetime(2003, 10, 1, 0, 0)

time_array_one_row = np.linspace(0, 8760, 8761)

time_array_datetime = [
    start_datetime + timedelta(hours=float(i))
    for i in time_array_one_row
]

# --------------------------------------------------
# 3. Load model outputs
# --------------------------------------------------
no_cap = np.load("little_washita_time72hr_sim_typeno-cap_storage_cycle5.npz")
subsurf_storage_no_cap = no_cap["arr_0"]
surf_storage_no_cap    = no_cap["arr_1"]
et_no_cap              = no_cap["arr_2"] #m3/hr
overland_no_cap        = no_cap["arr_3"]

cap = np.load("little_washita_time72hr_sim_typecap_storage.npz")
subsurf_storage_cap = cap["arr_0"]
surf_storage_cap    = cap["arr_1"]
et_cap              = cap["arr_2"] #m3/hr
overland_cap        = cap["arr_3"]

# --------------------------------------------------
# 4. Match lengths
# --------------------------------------------------
n = min(
    len(time_array_datetime),
    len(et_cap),
    len(et_no_cap),
    len(subsurf_storage_cap),
    len(subsurf_storage_no_cap),
    len(surf_storage_cap),
    len(surf_storage_no_cap)
)

time_array_datetime = time_array_datetime[:n]
et_cap = et_cap[:n]
et_no_cap = et_no_cap[:n]
subsurf_storage_cap = subsurf_storage_cap[:n]
subsurf_storage_no_cap = subsurf_storage_no_cap[:n]
surf_storage_cap = surf_storage_cap[:n]
surf_storage_no_cap = surf_storage_no_cap[:n]

# --------------------------------------------------
# 5. Convert model hourly ET to daily ET depth
# --------------------------------------------------
model_df = pd.DataFrame({
    "time": pd.to_datetime(time_array_datetime),
    "ET_cap_m3_hr": et_cap,
    "ET_no_cap_m3_hr": et_no_cap
})

# If model ET is m3/hr sampled hourly:
# summing hourly values over a day gives m3/day.
model_daily = (
    model_df
    .set_index("time")
    .resample("D")
    .sum()
    .reset_index()
)

model_daily["ET_cap_m3_day"] = model_daily["ET_cap_m3_hr"]
model_daily["ET_no_cap_m3_day"] = model_daily["ET_no_cap_m3_hr"]

model_daily["ET_cap_m_day"] = (
    model_daily["ET_cap_m3_day"] / (model_area_m2)
)

model_daily["ET_no_cap_m_day"] = (
    model_daily["ET_no_cap_m3_day"] / (model_area_m2)
)

# --------------------------------------------------
# 6. Load SSEBop daily mean ET
# --------------------------------------------------
# This file should contain regional daily mean ET in mm/day.
ssebop = pd.read_csv("ssebop_eta_conus2_bbox_daily_mean.csv")

ssebop["time"] = pd.to_datetime(ssebop["time"])

ssebop = (
    ssebop
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

# Convert SSEBop mm/day to m/day
ssebop["SSEBop_ET_m_day"] = ssebop["ET_mm_day"] / 1000.0

# --------------------------------------------------
# 7. Merge model and SSEBop by date
# --------------------------------------------------
compare = pd.merge(
    model_daily[
        [
            "time",
            "ET_cap_m3_day",
            "ET_no_cap_m3_day",
            "ET_cap_m_day",
            "ET_no_cap_m_day"
        ]
    ],
    ssebop[
        [
            "time",
            "ET_mm_day",
            "SSEBop_ET_m_day"
        ]
    ],
    on="time",
    how="inner"
)

compare = (
    compare
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

compare.to_csv(
    "model_vs_ssebop_daily_ET_mean_depth.csv",
    index=False
)

print("\nMatched comparison table:")
print(compare.head())
print("\nNumber of matched unique days:", len(compare))

# --------------------------------------------------
# 8. Statistics function
# --------------------------------------------------
def comparison_stats(model, obs):
    model = np.asarray(model, dtype=float)
    obs = np.asarray(obs, dtype=float)

    valid = np.isfinite(model) & np.isfinite(obs)

    model = model[valid]
    obs = obs[valid]

    bias = np.mean(model - obs)
    rmse = np.sqrt(np.mean((model - obs) ** 2))
    mae = np.mean(np.abs(model - obs))

    if len(model) > 1:
        corr = np.corrcoef(model, obs)[0, 1]
    else:
        corr = np.nan

    pbias = 100.0 * np.sum(model - obs) / np.sum(obs)

    return bias, rmse, mae, corr, pbias

# --------------------------------------------------
# 9. Statistics in m/day
# --------------------------------------------------
bias_cap, rmse_cap, mae_cap, r_cap, pbias_cap = comparison_stats(
    compare["ET_cap_m_day"],
    compare["SSEBop_ET_m_day"]
)

bias_no_cap, rmse_no_cap, mae_no_cap, r_no_cap, pbias_no_cap = comparison_stats(
    compare["ET_no_cap_m_day"],
    compare["SSEBop_ET_m_day"]
)

stats_df = pd.DataFrame({
    "case": ["Capillarity", "No capillarity"],
    "bias_m_day": [bias_cap, bias_no_cap],
    "rmse_m_day": [rmse_cap, rmse_no_cap],
    "mae_m_day": [mae_cap, mae_no_cap],
    "correlation": [r_cap, r_no_cap],
    "percent_bias": [pbias_cap, pbias_no_cap],
    "total_model_ET_m": [
        compare["ET_cap_m_day"].sum(),
        compare["ET_no_cap_m_day"].sum()
    ],
    "total_SSEBop_ET_m": [
        compare["SSEBop_ET_m_day"].sum(),
        compare["SSEBop_ET_m_day"].sum()
    ]
})

stats_df.to_csv(
    "model_vs_ssebop_daily_ET_mean_depth_stats.csv",
    index=False
)

print("\nComparison statistics using daily mean ET depth:")
print(stats_df)







# --------------------------------------------------
# 10. Plot daily ET depth comparison
# --------------------------------------------------
plt.figure(figsize=(14, 6), dpi=150)

plt.plot(
    compare["time"],
    compare["ET_cap_m_day"],
    linewidth=3,
    color=blue,
    linestyle="-",
    label="Model capillarity"
)

plt.plot(
    compare["time"],
    compare["ET_no_cap_m_day"],
    linewidth=3,
    color=red,
    linestyle="--",
    label="Model no capillarity"
)

plt.plot(
    compare["time"],
    compare["SSEBop_ET_m_day"],
    linewidth=3,
    color="k",
    linestyle="--",
    label="SSEBop"
)

plt.ylabel("Daily mean ET [m/day]")
plt.xlabel("Calendar date")
plt.legend(loc="best")
plt.xticks(rotation=30)
plt.tight_layout()

plt.savefig(
    "model_vs_ssebop_daily_ET_mean_depth.png",
    bbox_inches="tight",
    dpi=600
)

plt.savefig(
    "model_vs_ssebop_daily_ET_mean_depth.pdf",
    bbox_inches="tight",
    dpi=600
)

plt.show()

# --------------------------------------------------
# 11. Scatter comparison
# --------------------------------------------------
plt.figure(figsize=(7, 7), dpi=150)

plt.scatter(
    compare["SSEBop_ET_m_day"],
    compare["ET_cap_m_day"],
    s=40,
    color=blue,
    label="Capillarity"
)

plt.scatter(
    compare["SSEBop_ET_m_day"],
    compare["ET_no_cap_m_day"],
    s=40,
    color=red,
    label="No capillarity"
)

min_val = min(
    compare["SSEBop_ET_m_day"].min(),
    compare["ET_cap_m_day"].min(),
    compare["ET_no_cap_m_day"].min()
)

max_val = max(
    compare["SSEBop_ET_m_day"].max(),
    compare["ET_cap_m_day"].max(),
    compare["ET_no_cap_m_day"].max()
)

plt.plot(
    [min_val, max_val],
    [min_val, max_val],
    color=gray,
    linewidth=2,
    linestyle=":"
)

plt.xlabel("SSEBop ET [m/day]")
plt.ylabel("Model ET [m/day]")
plt.legend(loc="best")
plt.tight_layout()

plt.savefig(
    "model_vs_ssebop_daily_ET_mean_depth_scatter.png",
    bbox_inches="tight",
    dpi=600
)

plt.savefig(
    "model_vs_ssebop_daily_ET_mean_depth_scatter.pdf",
    bbox_inches="tight",
    dpi=600
)

plt.show()

# --------------------------------------------------
# 12. Water balance plot with daily mean ET panel
# --------------------------------------------------
fig, axs = plt.subplots(
    3,
    sharex=True,
    figsize=(12, 10),
    dpi=100
)

axs[0].plot(
    time_array_datetime,
    subsurf_storage_cap,
    linewidth=3,
    color=blue,
    linestyle="-",
    label="Capillarity"
)

axs[0].plot(
    time_array_datetime,
    subsurf_storage_no_cap,
    linewidth=3,
    color=red,
    linestyle="--",
    label="No capillarity"
)

axs[0].legend(loc="best")
axs[0].set_ylabel("Ss [m$^3$]")

axs[1].plot(
    time_array_datetime,
    surf_storage_cap,
    linewidth=3,
    color=blue,
    linestyle="-"
)

axs[1].plot(
    time_array_datetime,
    surf_storage_no_cap,
    linewidth=3,
    color=red,
    linestyle="--"
)

axs[1].set_ylabel("S [m$^3$]")


axs[2].plot(
    compare["time"],
    compare["ET_cap_m_day"]*model_area_m2,
    linewidth=3,
    color=blue,
    linestyle="-",
)

axs[2].plot(
    compare["time"],
    compare["ET_no_cap_m_day"]*model_area_m2,
    linewidth=3,
    color=red,
    linestyle="--",
)

axs[2].plot(
    compare["time"],
    compare["SSEBop_ET_m_day"]*model_area_m2,
    linewidth=3,
    color="k",
    linestyle="--",
    label="SSEBop"
)

axs[2].set_ylabel("Daily ET [m$^3$/day]")
axs[2].legend(loc="best")

plt.xlabel("Calendar date")

plt.subplots_adjust(
    left=0.1,
    bottom=0.1,
    right=0.9,
    top=0.85,
    wspace=0.0,
    hspace=0.0
)

axs[1].ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(0, 0)
)

axs[2].ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(0, 0)
)

plt.xticks(rotation=30)

plt.savefig(
    "little_washita_storage_ET_mean_depth_SSEBop.pdf",
    bbox_inches="tight",
    dpi=600
)

plt.savefig(
    "little_washita_storage_ET_mean_depth_SSEBop.png",
    bbox_inches="tight",
    dpi=600
)

plt.show()

# --------------------------------------------------
# 13. Print totals
# --------------------------------------------------
print("\nAll model totals:")

print(
    "With capillarity:",
    "Surface storage change",
    surf_storage_cap[-1] - surf_storage_cap[0],
    "Subsurface storage change",
    subsurf_storage_cap[-1] - subsurf_storage_cap[0],
    "ET total volume [m3]",
    compare["ET_cap_m3_day"].sum(),
    "ET total depth [m]",
    compare["ET_cap_m_day"].sum()
)

print(
    "Without capillarity:",
    "Surface storage change",
    surf_storage_no_cap[-1] - surf_storage_no_cap[0],
    "Subsurface storage change",
    subsurf_storage_no_cap[-1] - subsurf_storage_no_cap[0],
    "ET total volume [m3]",
    compare["ET_no_cap_m3_day"].sum(),
    "ET total depth [m]",
    compare["ET_no_cap_m_day"].sum()
)

print(
    "SSEBop:",
    "ET total depth [m]",
    compare["SSEBop_ET_m_day"].sum()
)