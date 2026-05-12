#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 07:49:47 2026

@author: ms6985
"""

import numpy as np
import matplotlib.pyplot as plt


plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'font.family': 'Serif'})

#Colors
brown  = [181/255 , 101/255, 29/255]
red    = [190/255 ,30/255 ,45/255 ]
blue   = [ 30/255 ,144/255 , 255/255 ]
green  = [  0/255 , 166/255 ,  81/255]
orange = [247/255 , 148/255 ,  30/255]
purple = [102/255 ,  45/255 , 145/255]
brown  = [155/255 ,  118/255 ,  83/255]
tan    = [199/255 , 178/255 , 153/255]
gray   = [100/255 , 100/255 , 100/255]

from datetime import datetime, timedelta
d = datetime(2003, 10, 1, 0, 0)

#print(d + timedelta(hours = time_array_plot[0,2]))
time_array_one_row = np.linspace(0,8760,8761)#time_array_plot[0,:]

time_array_datetime = ([d + timedelta(hours = i) for i in time_array_one_row])

#Difference plots  subsurface_storage,surface_storage,et,overland_flow
no_cap = np.load(f'little_washita_time72hr_sim_typeno-cap_storage_cycle5.npz'); subsurf_storage_no_cap = no_cap['arr_0']; surf_storage_no_cap = no_cap['arr_1']; et_no_cap = no_cap['arr_2']; overland_no_cap = no_cap['arr_3'];  
cap = np.load(f'little_washita_time72hr_sim_typecap_storage.npz'); subsurf_storage_cap = cap['arr_0']; surf_storage_cap = cap['arr_1']; et_cap = cap['arr_2']; overland_cap = cap['arr_3'];  

#plotting components of water balance
fig, axs = plt.subplots(3, sharex=True, figsize=(12,10) , dpi=100)
#fig.suptitle('Little Washita Water Balance')
axs[0].plot(time_array_datetime,subsurf_storage_cap, linewidth=3, color = blue,linestyle='-',label='Capillarity')
#ax2=axs[0].twinx()
axs[0].plot(time_array_datetime,subsurf_storage_no_cap, linewidth=3, color = red,linestyle='--',label='No capillarity')
axs[0].legend(loc='best')
axs[0].set_ylabel("Ss [m$^3$]")
axs[1].plot(time_array_datetime,surf_storage_cap, linewidth=3, color = blue,linestyle='-')
axs[1].plot(time_array_datetime,surf_storage_no_cap, linewidth=3, color = red,linestyle='--')
axs[1].set_ylabel("S [m$^3$]")
axs[2].plot(time_array_datetime,et_cap, linewidth=3, color = blue,linestyle='-')
axs[2].plot(time_array_datetime,et_no_cap, linewidth=3, color = red,linestyle='--')
axs[2].set_ylabel("ET [m$^3$/hr]")
plt.xlabel('Calendar date')
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.85, 
                    wspace=0.0, 
                    hspace=0.0)
axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xticks(rotation = 30)
plt.savefig(f'little_washita_time_year_sim_type_combined_no_Q_spun_up_abs.pdf',bbox_inches='tight', dpi = 600) 
plt.savefig(f'little_washita_time_year_sim_type_combined_no_Q_spun_up_abs.png',bbox_inches='tight', dpi = 600) 

print('All solutions with total in m3')
print('With capillarity: Surface storage',surf_storage_cap[-1]-subsurf_storage_cap[0],'Change in subsurface storage',subsurf_storage_cap[-1]-subsurf_storage_cap[0], 'ET',np.sum(et_cap))
print('Without capillarity: Surface storage',surf_storage_no_cap[-1]-subsurf_storage_no_cap[0],'Change in subsurface storage',subsurf_storage_no_cap[-1]-subsurf_storage_no_cap[0], 'ET',np.sum(et_no_cap))
print('Percentage: Surface storage',(surf_storage_no_cap[-1]-subsurf_storage_no_cap[0]-surf_storage_cap[-1]-subsurf_storage_cap[0])/(surf_storage_cap[-1]-subsurf_storage_cap[0])*100,'Change in subsurface storage',(subsurf_storage_no_cap[-1]-subsurf_storage_no_cap[0]-(subsurf_storage_cap[-1]-subsurf_storage_cap[0]))/(subsurf_storage_cap[-1]-subsurf_storage_cap[0])*100, 'ET',(np.sum(et_no_cap)-np.sum(et_cap))/(np.sum(et_cap))*100)









# compare_model_vs_ssebop_daily_et.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# --------------------------------------------------
# Plot settings
# --------------------------------------------------
plt.rcParams.update({"font.size": 20})
plt.rcParams.update({"font.family": "Serif"})

brown  = [181/255, 101/255, 29/255]
red    = [190/255,  30/255, 45/255]
blue   = [ 30/255, 144/255, 255/255]
green  = [  0/255, 166/255, 81/255]
orange = [247/255, 148/255, 30/255]
purple = [102/255,  45/255, 145/255]
tan    = [199/255, 178/255, 153/255]
gray   = [100/255, 100/255, 100/255]

# --------------------------------------------------
# 1. Time axis for model
# --------------------------------------------------
start_datetime = datetime(2003, 10, 1, 0, 0)

# 8760 hours + initial time = 8761 timestamps
time_array_one_row = np.linspace(0, 8760, 8761)

time_array_datetime = [
    start_datetime + timedelta(hours=float(i))
    for i in time_array_one_row
]

# --------------------------------------------------
# 2. Load model outputs
# --------------------------------------------------
no_cap = np.load("little_washita_time72hr_sim_typeno-cap_storage_cycle5.npz")
subsurf_storage_no_cap = no_cap["arr_0"]
surf_storage_no_cap    = no_cap["arr_1"]
et_no_cap              = no_cap["arr_2"]
overland_no_cap        = no_cap["arr_3"]

cap = np.load("little_washita_time72hr_sim_typecap_storage.npz")
subsurf_storage_cap = cap["arr_0"]
surf_storage_cap    = cap["arr_1"]
et_cap              = cap["arr_2"]
overland_cap        = cap["arr_3"]

# --------------------------------------------------
# 3. Make sure model arrays and time have same length
# --------------------------------------------------
n = min(len(time_array_datetime), len(et_cap), len(et_no_cap))

time_array_datetime = time_array_datetime[:n]
et_cap = et_cap[:n]
et_no_cap = et_no_cap[:n]

# --------------------------------------------------
# 4. Convert hourly model ET to daily ET
#    Model ET assumed to be m3/hr.
#    Daily sum gives m3/day.
# --------------------------------------------------
model_df = pd.DataFrame({
    "time": pd.to_datetime(time_array_datetime),
    "ET_cap_m3_hr": et_cap,
    "ET_no_cap_m3_hr": et_no_cap
})

model_df = (
    model_df
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

model_daily = (
    model_df
    .set_index("time")
    .resample("D")
    .sum()
    .reset_index()
)

model_daily["ET_cap_m3_day"] = model_daily["ET_cap_m3_hr"]*24
model_daily["ET_no_cap_m3_day"] = model_daily["ET_no_cap_m3_hr"]*24

# --------------------------------------------------
# 5. Load SSEBop daily summed ET
# --------------------------------------------------
ssebop = pd.read_csv("ssebop_eta_conus2_bbox_daily_sum.csv")

ssebop["time"] = pd.to_datetime(ssebop["time"])

ssebop = (
    ssebop
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

# --------------------------------------------------
# 6. Convert SSEBop to m3/day
# --------------------------------------------------
# Your SSEBop CSV has ET_mm_day = sum of mm/day over pixels.
# For 1 km x 1 km cells:
#cell_area_m2 = 1000.0 * 1000.0


# --------------------------------------------------
# 7. Merge model and SSEBop by unique date
# --------------------------------------------------
compare = pd.merge(
    model_daily[["time", "ET_cap_m3_day", "ET_no_cap_m3_day"]],
    ssebop[["time", "ET_m3_day"]],
    on="time",
    how="inner"
)

compare = (
    compare
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

compare.to_csv("model_vs_ssebop_daily_ET.csv", index=False)

print("\nMatched comparison table:")
print(compare.head())
print("\nNumber of matched unique days:", len(compare))

# --------------------------------------------------
# 8. Plot comparison
# --------------------------------------------------
plt.figure(figsize=(14, 6), dpi=150)

plt.plot(
    compare["time"],
    compare["ET_cap_m3_day"],
    linewidth=3,
    color=blue,
    linestyle="-",
    label="Model capillarity"
)

plt.plot(
    compare["time"],
    compare["ET_no_cap_m3_day"],
    linewidth=3,
    color=red,
    linestyle="--",
    label="Model no capillarity"
)

plt.plot(
    compare["time"],
    compare["ET_m3_day"],
    linewidth=3,
    color=green,
    linestyle="-.",
    label="SSEBop"
)

plt.ylabel("Daily ET [m$^3$/day]")
plt.xlabel("Calendar date")
plt.legend(loc="best")
plt.xticks(rotation=30)
plt.tight_layout()

plt.savefig("model_vs_ssebop_daily_ET.png", bbox_inches="tight", dpi=600)
plt.savefig("model_vs_ssebop_daily_ET.pdf", bbox_inches="tight", dpi=600)
plt.show()

# --------------------------------------------------
# 9. Scatter plots
# --------------------------------------------------
plt.figure(figsize=(7, 7), dpi=150)

plt.scatter(
    compare["ET_m3_day"],
    compare["ET_cap_m3_day"],
    s=40,
    color=blue,
    label="Capillarity"
)

plt.scatter(
    compare["ET_m3_day"],
    compare["ET_no_cap_m3_day"],
    s=40,
    color=red,
    label="No capillarity"
)

min_val = min(
    compare["ET_m3_day"].min(),
    compare["ET_cap_m3_day"].min(),
    compare["ET_no_cap_m3_day"].min()
)

max_val = max(
    compare["ET_m3_day"].max(),
    compare["ET_cap_m3_day"].max(),
    compare["ET_no_cap_m3_day"].max()
)

plt.plot(
    [min_val, max_val],
    [min_val, max_val],
    color=gray,
    linewidth=2,
    linestyle=":"
)

plt.xlabel("SSEBop ET [m$^3$/day]")
plt.ylabel("Model ET [m$^3$/day]")
plt.legend(loc="best")
plt.tight_layout()

plt.savefig("model_vs_ssebop_daily_ET_scatter.png", bbox_inches="tight", dpi=600)
plt.savefig("model_vs_ssebop_daily_ET_scatter.pdf", bbox_inches="tight", dpi=600)
plt.show()

# --------------------------------------------------
# 10. Summary statistics
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


bias_cap, rmse_cap, mae_cap, r_cap, pbias_cap = comparison_stats(
    compare["ET_cap_m3_day"],
    compare["ET_m3_day"]
)

bias_no_cap, rmse_no_cap, mae_no_cap, r_no_cap, pbias_no_cap = comparison_stats(
    compare["ET_no_cap_m3_day"],
    compare["ET_m3_day"]
)

stats_df = pd.DataFrame({
    "case": ["Capillarity", "No capillarity"],
    "bias_m3_day": [bias_cap, bias_no_cap],
    "rmse_m3_day": [rmse_cap, rmse_no_cap],
    "mae_m3_day": [mae_cap, mae_no_cap],
    "correlation": [r_cap, r_no_cap],
    "percent_bias": [pbias_cap, pbias_no_cap],
    "total_model_ET_m3": [
        compare["ET_cap_m3_day"].sum(),
        compare["ET_no_cap_m3_day"].sum()
    ],
    "total_SSEBop_ET_m3": [
        compare["ET_m3_day"].sum(),
        compare["ET_m3_day"].sum()
    ]
})

stats_df.to_csv("model_vs_ssebop_daily_ET_stats.csv", index=False)

print("\nComparison statistics:")
print(stats_df)

# --------------------------------------------------
# 11. Existing water-balance plot
# --------------------------------------------------
fig, axs = plt.subplots(3, sharex=True, figsize=(12, 10), dpi=100)

axs[0].plot(
    time_array_datetime,
    subsurf_storage_cap[:n],
    linewidth=3,
    color=blue,
    linestyle="-",
    label="Capillarity"
)

axs[0].plot(
    time_array_datetime,
    subsurf_storage_no_cap[:n],
    linewidth=3,
    color=red,
    linestyle="--",
    label="No capillarity"
)

axs[0].legend(loc="best")
axs[0].set_ylabel("Ss [m$^3$]")

axs[1].plot(
    time_array_datetime,
    surf_storage_cap[:n],
    linewidth=3,
    color=blue,
    linestyle="-"
)

axs[1].plot(
    time_array_datetime,
    surf_storage_no_cap[:n],
    linewidth=3,
    color=red,
    linestyle="--"
)

axs[1].set_ylabel("S [m$^3$]")

axs[2].plot(
    time_array_datetime,
    et_cap[:n],
    linewidth=3,
    color=blue,
    linestyle="-"
)

axs[2].plot(
    time_array_datetime,
    et_no_cap[:n],
    linewidth=3,
    color=red,
    linestyle="--"
)

axs[2].set_ylabel("ET [m$^3$/hr]")
plt.xlabel("Calendar date")

plt.subplots_adjust(
    left=0.1,
    bottom=0.1,
    right=0.9,
    top=0.85,
    wspace=0.0,
    hspace=0.0
)

axs[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.xticks(rotation=30)

plt.savefig(
    "little_washita_time_year_sim_type_combined_no_Q_spun_up_abs.pdf",
    bbox_inches="tight",
    dpi=600
)

plt.savefig(
    "little_washita_time_year_sim_type_combined_no_Q_spun_up_abs.png",
    bbox_inches="tight",
    dpi=600
)

plt.show()

# --------------------------------------------------
# 12. Print total model water-balance values
# --------------------------------------------------
print("\nAll solutions with total in m3")

print(
    "With capillarity:",
    "Surface storage",
    surf_storage_cap[-1] - surf_storage_cap[0],
    "Change in subsurface storage",
    subsurf_storage_cap[-1] - subsurf_storage_cap[0],
    "ET",
    np.sum(et_cap)
)

print(
    "Without capillarity:",
    "Surface storage",
    surf_storage_no_cap[-1] - surf_storage_no_cap[0],
    "Change in subsurface storage",
    subsurf_storage_no_cap[-1] - subsurf_storage_no_cap[0],
    "ET",
    np.sum(et_no_cap)
)



# --------------------------------------------------
# 13. Equivalent ET depth comparison: m/day
# --------------------------------------------------

model_width_m  = 64000.0
model_height_m = 32000.0
model_area_m2  = model_width_m * model_height_m

ssebop_width_m  = 63000.03
ssebop_height_m = 32674.32
ssebop_area_m2  = ssebop_width_m * ssebop_height_m

compare["ET_cap_m_day"] = compare["ET_cap_m3_day"] / model_area_m2
compare["ET_no_cap_m_day"] = compare["ET_no_cap_m3_day"] / model_area_m2
compare["SSEBop_ET_m_day"] = compare["ET_m3_day"] / ssebop_area_m2

compare.to_csv("model_vs_ssebop_daily_ET_m_day.csv", index=False)

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
    color=green,
    linestyle="-.",
    label="SSEBop"
)

plt.ylabel("Daily ET [m/day]")
plt.xlabel("Calendar date")
plt.legend(loc="best")
plt.xticks(rotation=30)
plt.tight_layout()

plt.savefig("model_vs_ssebop_daily_ET_m_day.png", bbox_inches="tight", dpi=600)
plt.savefig("model_vs_ssebop_daily_ET_m_day.pdf", bbox_inches="tight", dpi=600)
plt.show()

# --------------------------------------------------
# 14. Scatter comparison in m/day
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

plt.savefig("model_vs_ssebop_daily_ET_m_day_scatter.png", bbox_inches="tight", dpi=600)
plt.savefig("model_vs_ssebop_daily_ET_m_day_scatter.pdf", bbox_inches="tight", dpi=600)
plt.show()

# --------------------------------------------------
# 15. Statistics in m/day
# --------------------------------------------------
bias_cap_m, rmse_cap_m, mae_cap_m, r_cap_m, pbias_cap_m = comparison_stats(
    compare["ET_cap_m_day"],
    compare["SSEBop_ET_m_day"]
)

bias_no_cap_m, rmse_no_cap_m, mae_no_cap_m, r_no_cap_m, pbias_no_cap_m = comparison_stats(
    compare["ET_no_cap_m_day"],
    compare["SSEBop_ET_m_day"]
)

stats_m_day_df = pd.DataFrame({
    "case": ["Capillarity", "No capillarity"],
    "bias_m_day": [bias_cap_m, bias_no_cap_m],
    "rmse_m_day": [rmse_cap_m, rmse_no_cap_m],
    "mae_m_day": [mae_cap_m, mae_no_cap_m],
    "correlation": [r_cap_m, r_no_cap_m],
    "percent_bias": [pbias_cap_m, pbias_no_cap_m],
    "total_model_ET_m": [
        compare["ET_cap_m_day"].sum(),
        compare["ET_no_cap_m_day"].sum()
    ],
    "total_SSEBop_ET_m": [
        compare["SSEBop_ET_m_day"].sum(),
        compare["SSEBop_ET_m_day"].sum()
    ]
})

stats_m_day_df.to_csv("model_vs_ssebop_daily_ET_m_day_stats.csv", index=False)

print("\nComparison statistics in m/day:")
print(stats_m_day_df)









# --------------------------------------------------
# Prepare daily ET for panel 3
# --------------------------------------------------
import pandas as pd

# Hourly model ET -> dataframe
model_df = pd.DataFrame({
    "time": pd.to_datetime(time_array_datetime),
    "et_cap": et_cap,
    "et_no_cap": et_no_cap
})

# remove duplicate timestamps
model_df = (
    model_df
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

# Daily ET [m3/day]
model_daily = (
    model_df
    .set_index("time")
    .resample("D")
    .sum()
    .reset_index()
)

# SSEBop daily ET
ssebop = pd.read_csv("ssebop_eta_conus2_bbox_daily_sum.csv")

ssebop["time"] = pd.to_datetime(ssebop["time"])

ssebop = (
    ssebop
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

# --------------------------------------------------
# Plotting components of water balance
# --------------------------------------------------
fig, axs = plt.subplots(
    3,
    sharex=True,
    figsize=(12,10),
    dpi=100
)

# --------------------------------------------------
# Panel 1: subsurface storage
# --------------------------------------------------
axs[0].plot(
    time_array_datetime,
    subsurf_storage_cap,
    linewidth=3,
    color=blue,
    linestyle='-',
    label='Capillarity'
)

axs[0].plot(
    time_array_datetime,
    subsurf_storage_no_cap,
    linewidth=3,
    color=red,
    linestyle='--',
    label='No capillarity'
)

axs[0].legend(loc='best')
axs[0].set_ylabel("Ss [m$^3$]")

# --------------------------------------------------
# Panel 2: surface storage
# --------------------------------------------------
axs[1].plot(
    time_array_datetime,
    surf_storage_cap,
    linewidth=3,
    color=blue,
    linestyle='-'
)

axs[1].plot(
    time_array_datetime,
    surf_storage_no_cap,
    linewidth=3,
    color=red,
    linestyle='--'
)

axs[1].set_ylabel("S [m$^3$]")

# --------------------------------------------------
# Panel 3: DAILY ET [m3/day]
# --------------------------------------------------
axs[2].plot(
    model_daily["time"],
    model_daily["et_cap"],
    linewidth=3,
    color=blue,
    linestyle='-',
    label='Capillarity'
)

axs[2].plot(
    model_daily["time"],
    model_daily["et_no_cap"],
    linewidth=3,
    color=red,
    linestyle='--',
    label='No capillarity'
)

# SSEBop overlay
axs[2].plot(
    ssebop["time"],
    ssebop["ET_m3_day"],
    linewidth=3,
    color='k',
    linestyle='--',
    label='SSEBop'
)

axs[2].set_ylabel("ET [m$^3$/day]")
axs[2].legend(loc='best')

# --------------------------------------------------
# Formatting
# --------------------------------------------------
plt.xlabel('Calendar date')

plt.subplots_adjust(
    left=0.1,
    bottom=0.1,
    right=0.9,
    top=0.85,
    wspace=0.0,
    hspace=0.0
)

axs[1].ticklabel_format(
    axis='y',
    style='sci',
    scilimits=(0,0)
)

axs[2].ticklabel_format(
    axis='y',
    style='sci',
    scilimits=(0,0)
)

plt.xticks(rotation=30)

plt.savefig(
    'little_washita_time_year_sim_type_combined_ET_SSEBop.pdf',
    bbox_inches='tight',
    dpi=600
)

plt.savefig(
    'little_washita_time_year_sim_type_combined_ET_SSEBop.png',
    bbox_inches='tight',
    dpi=600
)

plt.show()







# --------------------------------------------------
# Prepare daily ET for panel 3
# --------------------------------------------------
import pandas as pd

# Hourly model ET -> dataframe
model_df = pd.DataFrame({
    "time": pd.to_datetime(time_array_datetime),
    "et_cap": et_cap,
    "et_no_cap": et_no_cap
})

# remove duplicate timestamps
model_df = (
    model_df
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

# Daily ET [m3/day]
model_daily = (
    model_df
    .set_index("time")
    .resample("D")
    .sum()
    .reset_index()
)

# SSEBop daily ET
ssebop = pd.read_csv("ssebop_eta_conus2_bbox_daily_sum.csv")

ssebop["time"] = pd.to_datetime(ssebop["time"])

ssebop = (
    ssebop
    .sort_values("time")
    .drop_duplicates(subset="time", keep="first")
)

# --------------------------------------------------
# Plotting components of water balance
# --------------------------------------------------
fig, axs = plt.subplots(
    3,
    sharex=True,
    figsize=(12,10),
    dpi=100
)

# --------------------------------------------------
# Panel 1: subsurface storage
# --------------------------------------------------
axs[0].plot(
    time_array_datetime,
    subsurf_storage_cap,
    linewidth=3,
    color=blue,
    linestyle='-',
    label='Capillarity'
)

axs[0].plot(
    time_array_datetime,
    subsurf_storage_no_cap,
    linewidth=3,
    color=red,
    linestyle='--',
    label='No capillarity'
)

axs[0].legend(loc='best')
axs[0].set_ylabel("Ss [m$^3$]")

# --------------------------------------------------
# Panel 2: surface storage
# --------------------------------------------------
axs[1].plot(
    time_array_datetime,
    surf_storage_cap,
    linewidth=3,
    color=blue,
    linestyle='-'
)

axs[1].plot(
    time_array_datetime,
    surf_storage_no_cap,
    linewidth=3,
    color=red,
    linestyle='--'
)

axs[1].set_ylabel("S [m$^3$]")

# --------------------------------------------------
# Panel 3: DAILY ET [m3/day]
# --------------------------------------------------
axs[2].plot(
    model_daily["time"],
    model_daily["et_cap"]*24,
    linewidth=3,
    color=blue,
    linestyle='-',
)

axs[2].plot(
    model_daily["time"],
    model_daily["et_no_cap"]*24,
    linewidth=3,
    color=red,
    linestyle='--',
)

# SSEBop overlay
axs[2].plot(
    ssebop["time"],
    ssebop["ET_m3_day"],
    linewidth=3,
    color='k',
    linestyle='--',
    label='SSEBop'
)

axs[2].set_ylabel("Daily ET [m$^3$/day]")
axs[2].legend(loc='best')

# --------------------------------------------------
# Formatting
# --------------------------------------------------
plt.xlabel('Calendar date')

plt.subplots_adjust(
    left=0.1,
    bottom=0.1,
    right=0.9,
    top=0.85,
    wspace=0.0,
    hspace=0.0
)

axs[1].ticklabel_format(
    axis='y',
    style='sci',
    scilimits=(0,0)
)

axs[2].ticklabel_format(
    axis='y',
    style='sci',
    scilimits=(0,0)
)

plt.xticks(rotation=30)

plt.savefig(
    'little_washita_time_year_sim_type_combined_ET_SSEBop.pdf',
    bbox_inches='tight',
    dpi=600
)

plt.savefig(
    'little_washita_time_year_sim_type_combined_ET_SSEBop.png',
    bbox_inches='tight',
    dpi=600
)

plt.show()


