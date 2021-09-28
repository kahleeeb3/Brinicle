#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on Mon Dec 7 2020

A discretized point model for simulating advection diffusion reaction

The reaction implemented here is a simple precipitation
The three species are named u, v, and w
The reactions are
1) u + v -> w

The rates and constants of the five reactions are
    k1 = k1 1/(M*s)

The equations then become:
du/dt = -k1*u*v 
dv/dt = -k1*u*v
dw/dt = +k2*u*v

TO DO:
    - add concentration dependence to diffusion
    - add w saturation to prevent overgrowth
    - parallelize the code for performance

@author: tompkinn
"""

# %% Imports

import math
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib import cm
import convert
# %% Functions

#  This is the advection reaction diffusion solver.
def ard_pde(up,vp,wp,du,dv,dw,xvel,yvel,dt,dx,dy):
# This numerically integrates using finite differences and the forward
# Euler method.
        
    # Set up zero flux bc's in x.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
    upy = np.column_stack((up[:,1],up,up[:,-2]))
    vpy = np.column_stack((vp[:,1],vp,vp[:,-2]))
    wpy = np.column_stack((wp[:,1],wp,wp[:,-2]))
    # Make boundaries zero flux in y.
    upxy = np.vstack((upy[1,:], upy, upy[-2,:]))
    vpxy = np.vstack((vpy[1,:], vpy, vpy[-2,:]))
    wpxy = np.vstack((wpy[1,:], wpy, wpy[-2,:]))
            
    # Perform finite differences. (Diffusion)
    # On axis terms for u
    # Calculate differences
    uxx=dt*(du/dx**2)*(upxy[:,2:]+upxy[:,0:-2]-2*upxy[:,1:-1])
    uxx=uxx[1:-1,:] # Remove extra rows
    # Calculate differences
    uyy=dt*(du/dy**2)*(upxy[2:,:]+upxy[0:-2,:]-2*upxy[1:-1,:])
    uyy=uyy[:,1:-1] # Remove extra columns
    
    # On axis terms for v
    # Calculate differences
    vxx=dt*(dv/dx**2)*(vpxy[:,2:]+vpxy[:,0:-2]-2*vpxy[:,1:-1])
    vxx=vxx[1:-1,:] # Remove extra rows
    # Calculate differences
    vyy=dt*(dv/dy**2)*(vpxy[2:,:]+vpxy[0:-2,:]-2*vpxy[1:-1,:])
    vyy=vyy[:,1:-1] # Remove extra columns
    
    # On axis terms for w
    # Calculate differences
    wxx=dt*(dw/dx**2)*(wpxy[:,2:]+wpxy[:,0:-2]-2*wpxy[:,1:-1])
    wxx=wxx[1:-1,:] # Remove extra rows
    # Calculate differences
    wyy=dt*(dw/dy**2)*(wpxy[2:,:]+wpxy[0:-2,:]-2*wpxy[1:-1,:])
    wyy=wyy[:,1:-1] # Remove extra columns
   
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15;  # Set ff=1 to turn the fudge-factor off.
    
    # Diagonal terms for u
    # Calculate differences
    uxy=dt*((du/(dx**2+dy**2))*(upxy[2:,2:]+upxy[0:-2,0:-2]-2*upxy[1:-1,1:-1]))/ff
    uyx=dt*((du/(dx**2+dy**2))*(upxy[0:-2,2:]+upxy[2:,0:-2]-2*upxy[1:-1,1:-1]))/ff
    
    # Diagonal terms for v
    # Calculate differences
    vxy=dt*((dv/(dx**2+dy**2))*(vpxy[2:,2:]+vpxy[0:-2,0:-2]-2*vpxy[1:-1,1:-1]))/ff
    vyx=dt*((dv/(dx**2+dy**2))*(vpxy[0:-2,2:]+vpxy[2:,0:-2]-2*vpxy[1:-1,1:-1]))/ff
    
    # Diagonal terms for w
    # Calculate differences
    wxy=dt*((dw/(dx**2+dy**2))*(wpxy[2:,2:]+wpxy[0:-2,0:-2]-2*wpxy[1:-1,1:-1]))/ff
    wyx=dt*((dw/(dx**2+dy**2))*(wpxy[0:-2,2:]+wpxy[2:,0:-2]-2*wpxy[1:-1,1:-1]))/ff
    
    # Caclulate the advection along the x axis.
    uxvel=dt*(vel_u*xvel/dx)*(upxy[0:-2,:]-upxy[1:-1,:])
    uxvel=uxvel[:,1:-1] # Remove extra rows
    vxvel=dt*(vel_v*xvel/dx)*(vpxy[0:-2,:]-vpxy[1:-1,:])
    vxvel=vxvel[:,1:-1] # Remove extra rows
    wxvel=dt*(vel_w*xvel/dx)*(wpxy[0:-2,:]-wpxy[1:-1,:])
    wxvel=wxvel[:,1:-1] # Remove extra rows
    
    # Caclulate the advection along the y axis.
    uyvel=dt*(vel_u*yvel/dy)*(upxy[:,0:-2]-upxy[:,1:-1])
    uyvel=uyvel[1:-1,:] # Remove extra rows
    vyvel=dt*(vel_v*yvel/dy)*(vpxy[:,0:-2]-vpxy[:,1:-1])
    vyvel=vyvel[1:-1,:] # Remove extra rows
    wyvel=dt*(vel_w*yvel/dy)*(wpxy[:,0:-2]-wpxy[:,1:-1])
    wyvel=wyvel[1:-1,:] # Remove extra rows
            
    # Combine diffusion with advection and kinetics using
    # the forward Euler algorithm.
    up=up+(uxx+uyy+uxy+uyx)+(uxvel+uyvel) + dt*(-k1*up*vp)
    vp=vp+(vxx+vyy+vxy+vyx)+(vxvel+vyvel) + dt*(-k1*up*vp)
    wp=wp+(wxx+wyy+wxy+wyx)+(wxvel+wyvel) + dt*(+k1*up*vp)
    
    # Here you apply any constant conditions
    # Holding the top left edge at a constant value
    tl_u = 1.0
    tl_v = 0.0
    # From fraction to end (top)
    up[0,int(np.round(5*ny/10)):] = tl_u
    vp[0,int(np.round(5*ny/10)):] = tl_v
    # Holding the bottom left edge at a constant value
    bl_u = 0.0
    bl_v = 1.0
    # From begining (bottom) to fraction
    up[0,0:int(np.round(5*ny/10))] = bl_u
    vp[0,0:int(np.round(5*ny/10))] = bl_v
        
    return [up, vp, wp]

# %% Model Parameters

# The time and space parameters for the model go here.
res = 0.5        # This sets the resolution for the simulation in mm/step.
Lx = 200         # This sets the x length of the simulation in mm.
Ly = 20          # This sets the y length of the simulation in mm.
tEnd = 2400      # This sets the duration of the simulation in s.
dt = 0.1         # The time step for the calculation in s.
dtWindow = 10    # This sets how often to update the plot in s.

# These are internal parameters for the simulation.
dx = res                           # The x resolution for the calculation.
dy = res                           # The y resolution for the calculation.
nx = math.floor(Lx/dx)             # The x dimension for storage.
ny = math.floor(Ly/dy)             # The y dimension for storage.
nt = math.floor(tEnd/dt)+1         # The t dimension for storage.
u = np.zeros((nx,ny,nt))           # Define the array for u data.
v = np.zeros((nx,ny,nt))           # Define the array for v data.
w = np.zeros((nx,ny,nt))           # Define the array for w data.

# %% Chemical Parameters

# The chemistry parameters of the model go here.
du = 0.002    # This sets the diffusion for u in mm^2/s.
dv = 0.002    # This sets the diffusion for v in mm^2/s.
dw = 0.000    # This sets the diffusion for v in mm^2/s.
xvel = 0.1    # This set the x advection velocity in mm/s.
yvel = 0.0    # This set the y advection velocity in mm/s.
vel_u = 1     # This adjusts the advection velocity for u.
vel_v = 1     # This adjusts the advection velocity for v.
vel_w = 0     # This adjusts the advection velocity for w.
k1 = 1.0      # Rate for the precipitation reaction in 1/(M*s).


# %% Initial Conditions

# This sets the initial conditions for the simulation. Initially every spot
# is set to the same values with some optional random variational
# of amplitude amp.
u0 = 0        # Arbitrary initial value
v0 = 0        # Arbitrary initial value
w0 = 0        # Arbitrary initial value
ampu = 0      # Set to zero to start everything the same
ampv = 0      # Set to zero to start everything the same
ampw = 0      # Set to zero to start everything the same
u[:,:,0] = u0 + ampu * np.random.rand(nx,ny) # Adds noise to IC
v[:,:,0] = v0 + ampv * np.random.rand(nx,ny) # Adds noise to IC
w[:,:,0] = w0 + ampw * np.random.rand(nx,ny) # Adds noise to IC
uC = 0.5      # The u value for the different spots.
vC = 0.5      # The v value for the different spots.
wC = 0.5      # The v value for the different spots.

# This is the middle spot of specified width
# spotwidth=5 # This is half width in steps
# spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
# spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
# spottop=int(np.round(ny/2))-spotwidth    # Determine the top edge
# spotbottom=int(np.round(ny/2))+spotwidth # Determine the bottom edge
# u[spotleft:spotright,spottop:spotbottom,0]=uC # Create the initial spot in u
# v[spotleft:spotright,spottop:spotbottom,0]=vC # Create the initial spot in v
# w[spotleft:spotright,spottop:spotbottom,0]=wC # Create the initial spot in w

# %% This section runs the simulation itself.
t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap < tEnd:    
       
    # From here on down is the nuts and bolts of the simulation.
    # Update u, v, and w for the next time output
    # using the pde function defined at the beginning of this file
    [u[:,:,t+1],v[:,:,t+1],w[:,:,t+1]]=ard_pde(u[:,:,t],v[:,:,t],w[:,:,t],du,dv,dw,xvel,yvel,dt,dx,dy)
    ups=u[1:4,1:4,t+1] # Samples for debugging
    vps=v[1:4,1:4,t+1] # Samples for debugging
    wps=w[1:4,1:4,t+1] # Samples for debugging
    t=t+1              # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    # This displays the step and time being plotted. 
    print('Step: {0} Time: {1:0.3f}s'.format(t,telap))
        
# %% Create the output files
    
files = []   # This is for the output files.
t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap <= tEnd:    

    # Make the plots
    x = np.arange(0,Lx,dx) # Create the x data for plotting
    y = np.arange(0,Ly,dy) # Create the y data for plotting
    Y, X = np.meshgrid(y,x) # Create the X,Y matrices for plotting
    Z1 = u[:,:,t] # Pull the Z data matrix for plotting
    Z2 = v[:,:,t] # Pull the Z data matrix for plotting
    Z3 = w[:,:,t] # Pull the Z data matrix for plotting
    Z4 = (u[:,:,t] + v[:,:,t]) # Pull the Z data matrix for plotting
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2) # Create the figure with subplots
    # Create the filled countour plot with colormap and manual levels
    cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm,levels=np.arange(0,1.1,0.01))
    fig.colorbar(cf1, ax=ax1) # Add the colorbar
    ax1.set_xlabel('x (mm)') # Label the x axis
    ax1.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration of u \n at {0:0.1f}s'.format(telap)
    ax1.set_title(utitle) # Title the plot
    ax1.set_aspect('equal') # Make the aspect ratio equal
    # Create the filled countour plot with colormap and manual levels
    cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm,levels=np.arange(0,1.1,0.01))
    fig.colorbar(cf2, ax=ax2) # Add the colorbar
    ax2.set_xlabel('x (mm)') # Label the x axis
    ax2.set_ylabel('y (mm)') # Label the y axis
    vtitle = 'Concentration of v \n at {0:0.1f}s'.format(telap)
    ax2.set_title(vtitle) # Title the plot
    ax2.set_aspect('equal') # Make the aspect ratio equal
    # Create the filled countour plot with colormap and manual levels
    cf3 = ax3.contourf(X,Y,Z3,cmap=cm.coolwarm,levels=np.arange(0,1.1,0.01))
    fig.colorbar(cf3, ax=ax3) # Add the colorbar
    ax3.set_xlabel('x (mm)') # Label the x axis
    ax3.set_ylabel('y (mm)') # Label the y axis
    wtitle = 'Concentration of w \n at {0:0.1f}s'.format(telap)
    ax3.set_title(wtitle) # Title the plot
    ax3.set_aspect('equal') # Make the aspect ratio equal
    # Create the filled countour plot with colormap and manual levels
    cf4 = ax4.contourf(X,Y,Z4,cmap=cm.coolwarm,levels=np.arange(0,2.2,0.01))
    fig.colorbar(cf4, ax=ax4) # Add the colorbar
    ax4.set_xlabel('x (mm)') # Label the x axis
    ax4.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration of u+v \n at {0:0.1f}s'.format(telap)
    ax4.set_title(utitle) # Title the plot
    ax4.set_aspect('equal') # Make the aspect ratio equal
    # plt.subplots_adjust(hspace=0.75,left=-0.05)
    plt.tight_layout()
    # plt.subplots_adjust(left=-0.3)

    # plt.show() # This shows the plots as the code is running
    fname = 'Precipitation/Precipitation_%06d.png' % t # Create the file name for each plot
    print('Saving frame', fname) # Print the status update
    fig.savefig(fname, dpi=300) # Save the image
    files.append(fname) # Update the filename
    plt.close(fig) # Close the image so it doesn't show while the code is running
    plt.clf() # Clear the figrue to save memory
    t=t+math.floor(dtWindow/dt)            # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    
#%% Export images to videos

convert.export("Precipitation")
