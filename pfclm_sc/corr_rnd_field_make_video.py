# Import the ParFlow package
from set_demo_defaults import *
from parflow import Run
import os
import shutil
from parflow.tools.fs import mkdir, cp, chdir, get_absolute_path, rm, exists


import matplotlib.animation as animation
#Make a video

import parflow.tools.hydrology as hydro
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from parflow.tools.fs import get_absolute_path, cp, rm, mkdir, exists
import parflow as pf

from datetime import datetime, timedelta
d = datetime(2020, 1, 1, 0, 0)

alpha_vG = 2 #Set alpha for vanGenuchten model 2 or 100 1/m
head_table = -2.5 #location of water table [m]

#print(d + timedelta(hours = time_array_plot[0,2]))
time_array_one_row = np.linspace(0,8760,8761)#time_array_plot[0,:]

t_time_datetime = ([d + timedelta(hours = i) for i in time_array_one_row])

time_array_datetime = np.tile(t_time_datetime, (20,1))

def surface_subsurface_storage_video(alpha_vG,dir_name,run_name):
# ============================================================
# MAKE A VIDEO OF SATURATION CROSS-SECTIONS USING contourf
# ============================================================

    run = Run.from_definition(f'{dir_name}/{run_name}.pfidb')
    data = run.data_accessor
    nt = len(data.times)
    nx = data.shape[2]
    ny = data.shape[1]
    nz = data.shape[0]
    dx = data.dx
    dy = data.dy
    dz = data.dz

    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)
    print(nx,dx,ny,dy,nz)
    
    #print(nt,nx,ny,nz,dx,dy,dz)
    
    porosity = data.computed_porosity 
    specific_storage = data.specific_storage 
    mannings = run.Mannings.Geom.domain.Value
    
    ## remove input filenames for TopoSlopes to force the data accessor to read the output slopes
    ## this fixes a windows issue
    run.TopoSlopesX.FileName = None
    run.TopoSlopesY.FileName = None
    
    slopex = data.slope_x 
    slopey = data.slope_y 
    mask   = data.mask
    
    # formatting the mask so that values outside the domain are NA and inside the domain are 1
    nanmask=mask.copy()
    nanmask[nanmask == 0] = 'NaN'
    nanmask[nanmask > 0] = 1
    print(slopex,slopey)
    plt.figure()
    plt.contourf(data.mask[:,0,:])
    plt.colorbar()

    files = glob(dir_name+f"/{run_name}.out.satur.*.pfb")
    sat = pf.read_pfb_sequence(files)
    
    # Prepare video writer
    writer = animation.FFMpegWriter(fps=10, bitrate=1800)
    
    video_name = f"saturation_video_alpha{alpha_vG}_head{head_table}m.mp4"
    print(f"Creating video: {video_name}")
    
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Set static labels
    ax.set_ylabel('z [m]')
    ax.set_xlabel('x [m]')
    
    # define color limits
    vmin = 0.1
    vmax = 1.0
    
        # ------------------------------------------------------------
    # Create edge-based coordinates for pcolormesh
    # (ParFlow fields are cell-centered; pcolormesh wants edges)
    # ------------------------------------------------------------
    x_edge = np.arange(0.0, (nx+1)*dx, dx)      # length nx+1
    z_edge = np.zeros(nz+1)
    z_edge[1:] = np.cumsum(dz)                  # length nz+1

    # Prepare video writer
    writer = animation.FFMpegWriter(fps=10, bitrate=1800)
    video_name = f"saturation_video_alpha{alpha_vG}.mp4"

    fig, ax = plt.subplots(figsize=(8, 6))
    print(f"Creating video: {video_name}")

    # Begin writing frames
    with writer.saving(fig, video_name, dpi=600):

        for icount in range(0,nt,50):

            ax.clear()

            ax.set_title(f"{t_time_datetime[icount]}")
            ax.set_ylabel("z [m]")
            ax.set_xlabel("x [m]")

            # ---- KEY LINE (your request) ----
            im = ax.pcolormesh(
                x_edge, z_edge,
                sat[icount, :, 0, :],
                vmin=0.1, vmax=1.0,
                cmap='Blues',
                shading='auto'
            )

            # Only add colorbar once
            if icount == 0:
                cbar = fig.colorbar(im, ax=ax)
                cbar.set_label("Saturation")

            writer.grab_frame()

            if icount % 100 == 0:
                print(f"Frame {icount}/{nt} complete...")
    
    print("Finished writing video.")
surface_subsurface_storage_video(alpha_vG,"/home/ms6985/ParFlow-project/pfclm_sc/output_water_balance","PFCLM_SC_runoff")

