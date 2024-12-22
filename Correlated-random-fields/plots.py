# Here we set up some functions for plotting outputs.
import parflow as pf
import numpy as np  # we will use numpy to work with the data
import matplotlib.pyplot as plt  # we will use matplotlib to plot the data
from parflow.tools.fs import get_absolute_path
from parflow.tools.io import write_pfb, read_pfb
from parflow import Run
import parflow.tools.hydrology as hydro
from glob import glob
import math
import os
from scipy.optimize import fsolve

from set_demo_defaults import *

def plot_domain(run_directory, variable, timestep=0):
    """Function to plot output from a ParFlow run"""

    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "domain_example.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz

    # Print a summary of the run data
    print(f"nx = {nx}, ny = {ny}, nz = {nz}, nt = {nt}")
    print(f"dx = {dx}, dy = {dy}, dz = {dz[0]}")

    # Load the data
    if variable == "porosity":
        data = read_pfb(get_absolute_path(f"domain_example.out.{variable}.pfb")).reshape(nz, nx)
    elif variable == "mannings":
        data = read_pfb(get_absolute_path(f"domain_example.out.mannings.pfb"))[0, :, :]
    else:
        data = read_pfb(get_absolute_path(f"domain_example.out.{variable}.{str(timestep).zfill(5)}.pfb")).reshape(nz, nx)
    
    # Set negative saturation values to NaN
    if variable == "satur":
        data[data < 0.0] = np.nan
    
    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)

    # Get limits for plotting
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    print(f"vmin: {vmin}, vmax: {vmax}")
    
    # Define labels for plots
    if variable == "satur":
        label = "Saturation [-]"
        title = "Saturation"
    elif variable == "press":
        label = "Pressure Head [m]"
        title = "Pressure Head"
    elif variable == "porosity":
        label = "Porosity"
        title = "Porosity"
    elif variable == "mannings":
        label = "Mannings"
        title = "Mannings"

    # Use pcolormesh to plot the data with the x and z coordinates with lines 
    # for the grid mesh from the ParFlow run grid
    fig, ax = plt.subplots()

    if variable == "mannings":
        im = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    else:
        im = ax.pcolormesh(x, z, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    plt.colorbar(im, ax=ax, label=label)
    
    # Include mesh lines
    if variable == "mannings":
        ax.hlines(y,x[0],x[-1],colors='white',linewidth=0.5)
        ax.vlines(x,y[0],y[-1],colors='white',linewidth=0.5)
    else:
        ax.hlines(z,x[0],x[-1],colors='white',linewidth=0.5)
        ax.vlines(x,z[0],z[-1],colors='white',linewidth=0.5)
    
    ax.set_xlabel('x [m]')
    if variable == "mannings":
        ax.set_ylabel('y [m]')
    else:
        ax.set_ylabel('z [m]')
    if variable in ["porosity", "mannings"]:
        ax.set_title(f"{title}")
    else:
        ax.set_title(f"{title} at t={timestep}")
    plt.show()


def plot_domain_corr_rnd(run_directory, variable, timestep=0):
    """Function to plot output from a ParFlow run"""

    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "harvey_flow.1.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    #nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz

    # Print a summary of the run data
    print(f"nx = {nx}, ny = {ny}, nz = {nz}")
    print(f"dx = {dx}, dy = {dy}, dz = {dz[0]}")

    # Load the data
    if variable == "porosity":
        data = read_pfb(get_absolute_path(f"harvey_flow.1.out.{variable}.pfb")).reshape(nz, nx)
    elif variable == "x-permeability":
        data = read_pfb(get_absolute_path(f"harvey_flow.1.out.perm_x.pfb")).reshape(nz, nx)
    elif variable == "z-permeability":
        data = read_pfb(get_absolute_path(f"harvey_flow.1.out.perm_z.pfb")).reshape(nz, nx)
    elif variable == "mannings":
        data = read_pfb(get_absolute_path(f"harvey_flow.1.out.mannings.pfb"))[0, :, :]
    else:
        data = read_pfb(get_absolute_path(f"harvey_flow.1.out.{variable}.{str(timestep).zfill(5)}.pfb")).reshape(nz, nx)
    
    # Set negative saturation values to NaN
    if variable == "satur":
        data[data < 0.0] = np.nan
    
    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)

    # Get limits for plotting
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    print(f"vmin: {vmin}, vmax: {vmax}")
    
    # Define labels for plots
    if variable == "satur":
        label = "Saturation [-]"
        title = "Saturation"
    elif variable == "press":
        label = "Pressure Head [m]"
        title = "Pressure Head"
    elif variable == "porosity":
        label = "Porosity"
        title = "Porosity"
    elif variable == "x-permeability":
        label = "X-Permeability"
        title = "X-Permeability"
    elif variable == "z-permeability":
        label = "Z-Permeability"
        title = "Z-Permeability"
    elif variable == "mannings":
        label = "Mannings"
        title = "Mannings"

    # Use pcolormesh to plot the data with the x and z coordinates with lines 
    # for the grid mesh from the ParFlow run grid
    fig, ax = plt.subplots()

    if variable == "mannings":
        im = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    else:
        im = ax.pcolormesh(x, z, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    plt.colorbar(im, ax=ax, label=label)
    
    # Include mesh lines
    if variable == "mannings":
        ax.hlines(y,x[0],x[-1],colors='white',linewidth=0.5)
        ax.vlines(x,y[0],y[-1],colors='white',linewidth=0.5)
    else:
        ax.hlines(z,x[0],x[-1],colors='white',linewidth=0.5)
        ax.vlines(x,z[0],z[-1],colors='white',linewidth=0.5)
    
    ax.set_xlabel('x [m]')
    if variable == "mannings":
        ax.set_ylabel('y [m]')
    else:
        ax.set_ylabel('z [m]')
    if variable in ["porosity", "mannings"]:
        ax.set_title(f"{title}")
    else:
        ax.set_title(f"{title} at t={timestep}")
    plt.show()



def plot_domain_mannings(run_directory):
    """Function to plot output from a ParFlow run: Mannings"""

    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "domain_example.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz

    # Print a summary of the run data
    print(f"nx = {nx}, ny = {ny}, nz = {nz}, nt = {nt}")
    print(f"dx = {dx}, dy = {dy}, dz = {dz[0]}")

    # Load the data
    data = read_pfb(get_absolute_path(f"domain_example.out.mannings.pfb"))[0, :, :]
    
    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)

    # Get limits for plotting
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    print(f"vmin: {vmin}, vmax: {vmax}")

    # Use pcolormesh to plot the data with the x and z coordinates with lines 
    # for the grid mesh from the ParFlow run grid
    fig, ax = plt.subplots()
    im = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    
    # Include mesh lines
    ax.hlines(z,x[0],x[-1],colors='white',linewidth=0.5)
    ax.vlines(x,z[0],z[-1],colors='white',linewidth=0.5)
    
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    plt.show()


def plot_vert_var(run_directory, variable, timestep=0):
    """Function to plot output from a ParFlow run"""

    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "domain_example.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz

    # Print a summary of the run data
    print(f"nx = {nx}, ny = {ny}, nz = {nz}, nt = {nt}")
    print(f"dx = {dx}, dy = {dy}, dz = {dz[0]}")

    # Load the data
    if variable == "porosity":
        data = read_pfb(get_absolute_path(f"domain_example.out.{variable}.pfb")).reshape(nz, nx)
    elif variable == "mannings":
        data = read_pfb(get_absolute_path(f"domain_example.out.mannings.pfb"))[0, :, :]
    else:
        data = read_pfb(get_absolute_path(f"domain_example.out.{variable}.{str(timestep).zfill(5)}.pfb")).reshape(nz, nx)
    
    # Set negative saturation values to NaN
    if variable == "satur":
        data[data < 0.0] = np.nan
    
    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)

    print(f"x = {x}, y = {y}, z = {z}")
    print(f"Shapes of : x = {x.shape}, y = {y.shape}, z = {z.shape}")

    # Get limits for plotting
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    print(f"vmin: {vmin}, vmax: {vmax}")
    
    # Define labels for plots
    if variable == "satur":
        label = "Saturation [-]"
        title = "Saturation"
    elif variable == "press":
        label = "Pressure Head [m]"
        title = "Pressure Head"
    elif variable == "porosity":
        label = "Porosity"
        title = "Porosity"
    elif variable == "mannings":
        label = "Mannings"
        title = "Mannings"

    # Use pcolormesh to plot the data with the x and z coordinates with lines 
    # for the grid mesh from the ParFlow run grid

    #im = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    fig = plt.figure(figsize=(8,8) , dpi=100)
    im = plt.plot(data,z[1:]-dz/2)#Plotting at cell centers
    plt.xlabel(f"{title}")
    plt.ylabel('z [m]')
    if variable in ["porosity", "mannings"]:
        plt.title(f"{title}")
    else:
        plt.title(f"{title} at t={timestep}")
    plt.show()



def plot_vert_var_combined(run_directory, variable, time_array,RelPerm_N,Saturation_N,alpha_vG,simulation_name):
    """Function to plot output from a ParFlow run"""

    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "domain_example.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz

    # Print a summary of the run data
    print(f"nx = {nx}, ny = {ny}, nz = {nz}, nt = {nt}")
    print(f"dx = {dx}, dy = {dy}, dz = {dz[0]}")


    
    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)

    print(f"x = {x}, y = {y}, z = {z}")
    print(f"Shapes of : x = {x.shape}, y = {y.shape}, z = {z.shape}")

    
    # Define labels for plots
    if variable == "satur":
        label = "Saturation [-]"
        title = "Saturation"
    elif variable == "press":
        label = "Pressure Head [m]"
        title = "Pressure Head"
    elif variable == "porosity":
        label = "Porosity"
        title = "Porosity"
    elif variable == "mannings":
        label = "Mannings"
        title = "Mannings"

    #Soil properties from input file
    Ks = 0.01465 #Saturated hydraulic conductivity = 0.01465 m/h
    phi = 0.25 #porosity (see Pfidb)
    #alpha_vG = 1.0
    s_s = 1.0; s_r = 0.2 #Saturation saturated and residual
    m = 1-1/RelPerm_N

    kr_vG= lambda h,n: (1-(np.abs(alpha_vG*h)**(n-1))/(1+np.abs(alpha_vG*h)**n)**(1-1/n))**2.0/((1 + np.abs(alpha_vG*h)**n)**((1-1/n)/2)) #kr, head in cms
    sw_vG= lambda h,n: (s_s - s_r)/((1 + np.abs(alpha_vG*h)**n)**(1-1/n))+s_r


    se_vG= lambda sw: ((sw - s_r)/(s_s - s_r))
    kr_sat_vG= lambda sw: se_vG(sw)**(1/2)*(1-(1-se_vG(sw)**(1/m))**m)**2
    
    lambda_vG= lambda sw: Ks/(phi*(s_s - s_r)) * (2*se_vG(sw)*((se_vG(sw))**(1/m)*(1 - (se_vG(sw))**(1/m))**m + ((se_vG(sw))**(1/m) - 1)*((1 - (se_vG(sw))**(1/m))**m - 1))*((1 - (se_vG(sw))**(1/m))**m - 1)/((se_vG(sw))**(1/m) - 1))
    
    #lambda_vG= lambda sw: Ks/(phi*(s_s - s_r)) *



    #lambda_vG = lambda sw: Ks*2*(se_vG(sw) - se_vG(sw)*(1 - se_vG(sw)**(1/m))**m) * (1-(1 - se_vG(sw)**(1/m))**m +(1 - se_vG(sw)**(1/m))**(m-1) * se_vG(sw)**(1/m)) *1/(phi*(s_s - s_r))

    #lambda_vG_numeric = lambda sw: Ks * ( kr_sat_vG(sw +1e-8) - kr_sat_vG(sw))/(phi*1e-8)

    dkrbydsw_vG_numeric = lambda sw: ( kr_sat_vG(sw +1e-8) - kr_sat_vG(sw))/(1e-8)  #numerical
    dkrbydsw_vG_analytic= lambda sw: 1/(s_s - s_r) * (0.5*(se_vG(sw))**0.5*((se_vG(sw))**(1/m) - 1)*((1 - (se_vG(sw))**(1/m))**m - 1) + 2*(se_vG(sw))**((0.5*m + 1)/m)*(1 - (se_vG(sw))**(1/m))**m)*((1 - (se_vG(sw))**(1/m))**m - 1)/((se_vG(sw))**1.0*((se_vG(sw))**(1/m) - 1))

    lambda_vG_numeric  = lambda sw: Ks/phi*dkrbydsw_vG_numeric(sw)

    lambda_vG_analytic = lambda sw: Ks/phi*dkrbydsw_vG_analytic(sw)

    fig = plt.figure(figsize=(5,7.5) , dpi=100)

    for timestep in time_array:
        # Load the data
        if variable == "porosity":
            data = read_pfb(get_absolute_path(f"domain_example.out.{variable}.pfb")).reshape(nz, nx)
        elif variable == "mannings":
            data = read_pfb(get_absolute_path(f"domain_example.out.mannings.pfb"))[0, :, :]
        else:
            data = read_pfb(get_absolute_path(f"domain_example.out.{variable}.{str(timestep).zfill(5)}.pfb")).reshape(nz, nx)
        
        # Set negative saturation values to NaN
        if variable == "satur":
            data[data < 0.0] = np.nan

        plt.plot(data,z[1:]-dz/2,'k-',label=f"t={timestep}",alpha=(time_array.index(timestep)+1)/(len(time_array)+1))


        if simulation_name == "wetting":
            #boundary_sat_init = 0.7
            Applied_flux = -0.005#Ks * kr_sat_vG(boundary_sat_init) #Applied flux = -0.001 m/h (downward direction)

            kr_vG_invert= lambda sw: Ks*kr_sat_vG(sw)+Applied_flux
            boundary_sat = fsolve(kr_vG_invert,0.9)  #calculate saturation using inversion

            Initial_head = 100 #m  initial head for the column (m)
            init_sat = sw_vG(Initial_head,Saturation_N) #initial saturation

            #boundary_sat = boundary_sat_init

            print('Initial saturation', init_sat)
            print('Boundary saturation', boundary_sat)

            print('Shock location',z[-1]+Applied_flux/(phi*(boundary_sat-init_sat))*timestep,'m') 
            plt.hlines(z[-1]+(Applied_flux)/(phi*(boundary_sat - init_sat))*timestep,init_sat,boundary_sat,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            plt.vlines(boundary_sat,z[-1],z[-1]+Applied_flux/(phi*(boundary_sat-init_sat))*timestep,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            plt.vlines(init_sat,z[0],z[-1]+Applied_flux/(phi*(boundary_sat-init_sat))*timestep,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))

        elif simulation_name == "unsaturated_column":
            s_s= 1.0 #Complete saturation
            init_sat = 0.7#sw_vG(Initial_head,Saturation_N)
            #Initial_head = 1 #m
            sw_vG= lambda h,n: (s_s - s_r)/((1 + abs(alpha_vG*h)**n)**(1-1/n))+s_r #sw, head in cms

            h_from_sat = lambda h: sw_vG(h,Saturation_N) - init_sat
            Initial_head = np.abs(fsolve(h_from_sat,1,xtol=1e-10))
            print(init_sat)
            sat_array = np.linspace(s_r,init_sat,10000)  #Saturation array
            print(timestep)
            lambda_vG_calc = lambda_vG(sat_array)  #Speed of rarefaction lambda_vG_numeric or lambda_vG
            lambda_vG_num_calc = lambda_vG_numeric(sat_array)  #Speed of rarefaction lambda_vG_numeric or lambda_vG

            print(np.linalg.norm(dkrbydsw_vG_numeric(sat_array)-dkrbydsw_vG_analytic(sat_array)))
            print(np.linalg.norm(lambda_vG_numeric(sat_array)-lambda_vG_analytic(sat_array)))
            #print(np.linalg.norm(lambda_vG_v2(sat_array)-lambda_vG(sat_array)))


            #plt.plot(sat_array, z[-1]-lambda_vG_calc*timestep,0,1,color=green,linestyle='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            #plt.plot(sat_array, z[-1]-lambda_vG_num_calc*timestep,0,1,color=red,linestyle='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            plt.plot(sat_array, z[-1]-lambda_vG_analytic(sat_array)*timestep,0,1,color=red,linestyle='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
   
            #print('Rising perch water table location',z[0]+(phi*(data[-1,0]-s_s))*timestep,'m') 
            plt.hlines(z[0]+(-kr_vG(Initial_head,RelPerm_N)*Ks)/(phi*(sw_vG(Initial_head,RelPerm_N)-s_s))*timestep, init_sat,s_s,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            plt.vlines(s_s,z[0],z[0]+(-kr_vG(Initial_head,RelPerm_N)*Ks)/(phi*(sw_vG(Initial_head,RelPerm_N)-s_s))*timestep,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            plt.vlines(init_sat,z[0]+(-kr_vG(Initial_head,RelPerm_N)*Ks)/(phi*(sw_vG(Initial_head,RelPerm_N)-s_s))*timestep,z[-1]-lambda_vG_numeric(init_sat)*timestep,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            
            #plt.vlines(data[-1,0],z[-1],z[-1]+Applied_flux/(0.25*(data[-1,0]-0.2))*timestep,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
            #plt.vlines(data[0,-1],z[0],z[-1]+Applied_flux/(0.25*(data[-1,0]-0.2))*timestep,color=red,linestyles='--',alpha=(time_array.index(timestep)+1)/(len(time_array)+1))
    plt.xlim([s_r-0.01,s_s+0.01])
    plt.xlabel(f"{title}")
    plt.ylabel('z [m]')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.legend(loc='best', frameon=False)
    plt.savefig(f'{simulation_name}_NRelPerm{RelPerm_N}_NSaturation{Saturation_N}_alphavG{alpha_vG}_s_w_vs_depth.pdf',bbox_inches='tight', dpi = 600)
    plt.show()



def plot_domain_mannings(run_directory):
    """Function to plot output from a ParFlow run: Mannings"""

    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "domain_example.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz

    # Print a summary of the run data
    print(f"nx = {nx}, ny = {ny}, nz = {nz}, nt = {nt}")
    print(f"dx = {dx}, dy = {dy}, dz = {dz[0]}")

    # Load the data
    data = read_pfb(get_absolute_path(f"domain_example.out.mannings.pfb"))[0, :, :]
    
    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)

    # Get limits for plotting
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    print(f"vmin: {vmin}, vmax: {vmax}")

    # Use pcolormesh to plot the data with the x and z coordinates with lines 
    # for the grid mesh from the ParFlow run grid
    fig, ax = plt.subplots()
    im = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    
    # Include mesh lines
    ax.hlines(z,x[0],x[-1],colors='white',linewidth=0.5)
    ax.vlines(x,z[0],z[-1],colors='white',linewidth=0.5)
    
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    plt.show()



def plot_domain_overland(run_directory, variable, timestep=0):
    """Function to plot output from a ParFlow run: Overland Flow module"""

    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "overland_tiltedV.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz

    # Print a summary of the run data
    print(f"nx = {nx}, ny = {ny}, nz = {nz}, nt = {nt}")
    print(f"dx = {dx}, dy = {dy}, dz = {dz[0]}")

    # Load the data
    data = read_pfb(get_absolute_path(f"overland_tiltedV.out.{variable}.{str(timestep).zfill(5)}.pfb"))[0, :, :]
    
    # Set negative saturation values to NaN
    if variable == "satur":
        data[data < 0.0] = np.nan
    
    # Set up x and z to match the shape of the ParFlow grid
    x = np.arange(0.0,(nx+1)*dx,dx)
    y = np.arange(0.0,(ny+1)*dy,dy)
    z = np.zeros(nz+1)
    z[1:] = np.cumsum(dz)

    # Get limits for plotting
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    print(f"vmin: {vmin}, vmax: {vmax}")
    
    # Define labels for plots
    if variable == "satur":
        label = "Saturation [-]"
        title = "Saturation"
    elif variable == "press":
        label = "Pressure Head [m]"
        title = "Pressure Head"

    # Use pcolormesh to plot the data with the x and z coordinates with lines 
    # for the grid mesh from the ParFlow run grid
    fig, ax = plt.subplots()
    im = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap='plasma_r')
    plt.colorbar(im, ax=ax, label=label)
    
    # Include mesh lines
    ax.hlines(z,x[0],x[-1],colors='white',linewidth=0.5)
    ax.vlines(x,z[0],z[-1],colors='white',linewidth=0.5)
    
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title(f"{title} at t={timestep}")
    plt.show()
    
def plot_subsurface_storage(run_directory, subset_slice=None):
    """Function to plot total subsurface storage over time based on a ParFlow run"""
    
    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "domain_example.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size n

    if subset_slice is None:
        subset_array = (slice(0, nz), slice(0, ny), slice(0, nx))
    else: 
        subset_array = (slice(subset_slice[0], subset_slice[1]),
                        slice(subset_slice[2], subset_slice[3]),
                        slice(subset_slice[4], subset_slice[5]))
        dz = dz[subset_slice[0]:subset_slice[1]]

    mask = data.mask[subset_array]
    porosity = data.computed_porosity[subset_array]
    specific_storage = data.specific_storage[subset_array]
    
    # Initialize empty array
    subsurface_storage = np.zeros(nt)
    
    press_files = sorted(glob(f'{run_directory}/domain_example.out.press*.pfb'))
    satur_files = sorted(glob(f'{run_directory}/domain_example.out.satur*.pfb'))
        
    # Iteratively calculate overland flow for whole domain
    for hour in range(0, nt): 

        pressure = pf.read_pfb(press_files[hour])[subset_array]
        saturation = pf.read_pfb(satur_files[hour])[subset_array]
            
        subsurface_storage[hour] = np.sum(hydro.calculate_subsurface_storage(porosity, pressure, saturation, specific_storage, dx, dy, dz, mask=mask),
        axis=(0, 1, 2))

    # Plot
    plt.plot(subsurface_storage)
    plt.xlabel("Time")
    plt.ylabel("Subsurface Storage")
    plt.title(f"Subsurface Storage over Time")
    plt.show()

def plot_streamflow(run_directory, grid_cell=None):
    """Function to plot a hydrograph for a single grid cell location based on a ParFlow run"""
    
    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, "overland_tiltedV.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz
    mask = data.mask
    slopex = (data.slope_x).squeeze()
    slopey = (data.slope_y).squeeze()

    # Initialize empty array
    olf = np.zeros((nt, ny, nx))
    
    press_files = sorted(glob(f'{run_directory}/overland_tiltedV.out.press*.pfb'))
    satur_files = sorted(glob(f'{run_directory}/overland_tiltedV.out.satur*.pfb'))
    mannings = (read_pfb(f"{run_directory}/overland_tiltedV.out.mannings.pfb")).squeeze()

    # Iteratively calculate overland flow for whole domain
    for hour in range(0, nt):        
        pressure = pf.read_pfb(press_files[hour])
        saturation = pf.read_pfb(satur_files[hour])
            
        olf[hour, ...] = hydro.calculate_overland_flow_grid(pressure, slopex, slopey, mannings, dx, dy, mask=mask)  

    # Extract a single location point for plotting. Default to halfway along x-axis and with y=0 unless otherwise noted
    if grid_cell == None:
        grid_cell = (math.floor(nx/2), 0)
        
    output_point = olf[:, grid_cell[0], grid_cell[1]]

    # Plot
    plt.plot(output_point)
    plt.xlabel("Time")
    plt.ylabel("Streamflow")
    plt.title(f"Streamflow at grid cell {grid_cell}")
    plt.show()

def plot_water_balance(run_directory, run_name):
    """Function to plot total subsurface storage over time based on a ParFlow run"""
    
    # Load the run from the file, this is the same as the run defined above
    run = Run.from_definition(os.path.join(run_directory, f"{run_name}.pfidb"))   

    data = run.data_accessor # get the data accessor, this makes it easier to access the data from the run
    nt = len(data.times)  # get the number of time steps
    nx = data.shape[2]    # get the number of cells in the x direction
    ny = data.shape[1]    # get the number of cells in the y direction
    nz = data.shape[0]    # get the number of cells in the z direction
    dx = data.dx          # get the cell size in the x direction
    dy = data.dy          # get the cell size in the y direction
    dz = data.dz          # get the cell size in the z direction, this is a 1D array of size nz
    mask = data.mask
    surface_mask = mask[-1, :, :].astype(int)
    porosity = data.computed_porosity
    specific_storage = data.specific_storage
    mannings = data.mannings
    slopex = data.slope_x 
    slopey = data.slope_y 

    surface_storage = np.zeros(nt)
    subsurface_storage = np.zeros(nt)
    et = np.zeros(nt)
    precip = np.zeros(nt)
    overland_out = np.zeros(nt)
    swe = np.zeros(nt)
    water_balance_ = np.zeros(nt)
    delta_storage = np.zeros(nt)


    for t in range(0, nt):
    
        data.time = t
        pressure = data.pressure 
        saturation = data.saturation

        subsurface_storage[t] = np.sum(hydro.calculate_subsurface_storage(porosity, pressure, saturation, specific_storage, dx,dy,dz, mask=mask),axis=(0,1,2))
        surface_storage[t] = np.sum(hydro.calculate_surface_storage(pressure, dx, dy, mask=mask), axis=(0,1))
        overland_out[t] = hydro.calculate_overland_flow(pressure, slopex, slopey, mannings, dx, dy, mask=mask)

        #generate array of rainfall from ParFlow keys
        # if (t <= run.Cycle.rainrec.r0.Length):
        #     precip[t]=-(run.Patch.z_upper.BCPressure.r0.Value)*(nx*dx)*(ny*dy)*run.TimeStep.Value  #  m/h over domain
        precip[t] = 0.01 * (nx*dx)*(ny*dy) * run.TimeStep.Value
        
        if t > 0: 
            data.forcing_time = t-1
            # precip[t] = np.sum(data.clm_forcing('APCP')[surface_mask==1])*(3600/1000)*(dx)*(dy)
            et[t]  = np.sum(data.clm_output('qflx_evap_tot')[surface_mask==1])*(3600/1000)*(dx)*(dy)
            swe[t]  = np.sum(data.clm_output('swe_out')[surface_mask==1])*(1/1000)*(dx)*(dy)

        if t > 1:
            delta_storage[t] = (subsurface_storage[t] + surface_storage[t] + swe[t]) - (subsurface_storage[t-1] + surface_storage[t-1] + swe[t-1])
            water_balance_[t] = -delta_storage[t] + precip[t] - et[t] - overland_out[t]
            
    # Plot
    plt.plot(water_balance_)
    plt.xlabel("Time")
    plt.ylabel("Water Balance")
    plt.title("Water Balance over Time")
    plt.show()

    plt.plot(delta_storage)
    plt.xlabel("Time")
    plt.ylabel("Change in Storage")
    plt.title("Change in Storage over Time")
    plt.show()