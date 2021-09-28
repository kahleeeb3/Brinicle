#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  25 2021

A reaction diffusion implementation of the Tyson-Fife model of BZ

@author: tompkinn
"""

# This program simulates diffusion

# %% Imports

import math
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib import cm
import convert

# %% Functions

# This is the steady state solver.
def SteadyState(x):
    # Returns u,v for no diffusion.
    up = x[0]
    vp = x[1]
    # These two lines implement the Brusselator.
    upp = c * up * (1 - up) - b * c * vp * (up - a) / (up + a)
    vpp = up - vp
    y = [upp,vpp]
    
    return y

# %% This is the reaction diffusion solver.
def pde(up,vp,du,dv,xvel,yvel,dt,dx,dy):
# This numerically integrates using finite differences and the forward
# Euler method.

    # Set up periodic bc's in y.
    # These pad the array with the last items from the other side,
    # Thus the last item wraps around to the other side.
    upy = np.column_stack((up[:,-1],up,up[:,1]))
    vpy = np.column_stack((vp[:,-1],vp,vp[:,1]))
    # Make boundaries periodic in x.
    upxy = np.vstack((upy[-1,:], upy, upy[1,:]))
    vpxy = np.vstack((vpy[-1,:], vpy, vpy[1,:]))
        
    # Set up zero flux bc's in x.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
#    upy = np.column_stack((up[:,1],up,up[:,-2]))
#    vpy = np.column_stack((vp[:,1],vp,vp[:,-2]))
    # Make boundaries zero flux in y.
#    upxy = np.vstack((upy[1,:], upy, upy[-2,:]))
#    vpxy = np.vstack((vpy[1,:], vpy, vpy[-2,:]))
            
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
    
    # Caclulate the advection along the x axis.
    uxvel=dt*(xvel/dx)*(upxy[0:-2,:]-upxy[1:-1,:])
    uxvel=uxvel[:,1:-1] # Remove extra rows
    vxvel=dt*(xvel/dx)*(vpxy[0:-2,:]-vpxy[1:-1,:])
    vxvel=vxvel[:,1:-1] # Remove extra rows
    
    # Caclulate the advection along the y axis.
    uyvel=dt*(yvel/dy)*(upxy[:,0:-2]-upxy[:,1:-1])
    uyvel=uyvel[1:-1,:] # Remove extra rows
    vyvel=dt*(yvel/dy)*(vpxy[:,0:-2]-vpxy[:,1:-1])
    vyvel=vyvel[1:-1,:] # Remove extra rows
            
    # Combine diffusion with advection and kinetics using
    # the forward Euler algorithm.
    up=up+(uxx+uyy+uxy+uyx)+(uxvel+uyvel) + dt*(c*up*(1-up)-b*c*vp*(up-a)/(up+a))
    vp=vp+(vxx+vyy+vxy+vyx)+(vxvel+vyvel) + dt*(up-vp)
    
    # Here you apply any constant conditions
    # Holding the top left edge at a constant value
    # tl_u = 1.0
    # tl_v = 0.0
    # up[0,int(np.round(ny/2)):] = tl_u
    # vp[0,int(np.round(ny/2)):] = tl_v
    # Holding the bottom left edge at a constant value
    # bl_u = 0.0
    # bl_v = 1.0
    # up[0,0:int(np.round(ny/2))] = bl_u
    # vp[0,0:int(np.round(ny/2))] = bl_v
        
    return [up, vp]

# %% Model Parameters

# The time and space parameters for the model go here.
res = 0.2        # This sets the resolution for the simulation in mm/step.
Lx = 20         # This sets the x length of the simulation in mm.
Ly = 10         # This sets the y length of the simulation in mm.
tEnd = 10       # This sets the duration of the simulation in s.
dt = 0.001        # The time step for the calculation in s.
dtWindow = 0.1   # This sets how often to update the plot in s.

# These are internal parameters for the simulation.
dx = res                           # The x resolution for the calculation.
dy = res                           # The y resolution for the calculation.
nx = math.floor(Lx/dx)             # The x dimension for storage.
ny = math.floor(Ly/dy)             # The y dimension for storage.
nt = math.floor(tEnd/dt)+1         # The t dimension for storage.
displayUandV = 1                   # Whether to plot u and v or just v.
u = np.zeros((nx,ny,nt))           # Define the array for u data.
v = np.zeros((nx,ny,nt))           # Define the array for v data.

# %% Chemical Parameters

"""
Try with Advection - Can you get patterns with flow?
"""

# The chemistry parameters of the model go here.
du = 0.002    # This sets the diffusion for u in mm^2/s.
dv = 1.0*du     # This sets the diffusion for v in mm^2/s.
xvel = 0.0    # This set the x advection velocity in mm/s.
yvel = 0.0    # This set the y advection velocity in mm/s.
a = 0.01     # This is parameter a from the model. a<<1
b = 1.49     # This is parameter b from the model. b~1
c = 100      # This is 1/eps from the model. eps<<1 thus c>>1

# %% Steady State Solution

# The code starts by finding the steady state values in the absence of
# diffusion. This is the fixed point for the system.
initGuess = [1,1]
uv0 = fsolve(SteadyState,initGuess)
u0 = uv0[0]   # Steady state value for u
v0 = uv0[1]   # Steady state value for v

# %% Initial Conditions

# This sets the initial conditions for the simulation. Initially every spot
# is set to the same values with some optional random variational
# of amplitude amp.
# u0 = 0.3        # Arbitrary initial value
# v0 = 0.3        # Arbitrary initial value
"""
Heterogeneity is also off - This model does not like it
have not played with it as much as the Brusselator
very small Heterogeneity - do you still get patterns?
How can you make Spirals?
"""
u_amp = u0/4      # Set to zero to start everything the same
v_amp = v0/4     # Set to zero to start everything the same
u[:,:,0] = u0*1 + u_amp * np.random.rand(nx,ny)
v[:,:,0] = v0*1 + v_amp * np.random.rand(nx,ny)
uC = u0/4          # The u value for the different spots.
vC = v0/4          # The v value for the different spots.

# This is the middle spot of specified width
spotwidth=5 # This is half width in steps
spotleft=int(np.round(nx/2))-spotwidth # Determine the left edge
spotright=int(np.round(nx/2))+spotwidth # Determine the right edge
spottop=int(np.round(ny/2))-spotwidth # Determine the top edge
spotbottom=int(np.round(ny/2))+spotwidth # Determine the bottom edge
u[spotleft:spotright,spottop:spotbottom,0]=uC # Create the initial spot in u
v[spotleft:spotright,spottop:spotbottom,0]=vC # Create the initial spot in y

# %% This section runs the simulation itself.
t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap < tEnd:    
       
    # From here on down is the nuts and bolts of the simulation.
    # Update u and v for the next time output point defined by nWindow
    # using the pde function defined at the beginning of this file
    [u[:,:,t+1],v[:,:,t+1]]=pde(u[:,:,t],v[:,:,t],du,dv,xvel,yvel,dt,dx,dy)
    ups=u[1:4,1:4,t+1] # Samples for debugging
    vps=v[1:4,1:4,t+1] # Samples for debugging
    t=t+1              # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    # This displays the step and time being plotted.    
    print('Step: {0} Time: {1}'.format(t,telap))
        
# %% Create the output files
    
files = []   # This is for the output files.
t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
plot_u_min = 0
plot_u_max = 1.0
plot_v_min = 0
plot_v_max = 0.3

while telap <= tEnd:    

    # Make the plots
    if displayUandV: # If the Boolean is set to plot both u and v
        x = np.arange(0,Lx,dx) # Create the x data for plotting
        y = np.arange(0,Ly,dy) # Create the y data for plotting
        Y, X = np.meshgrid(y,x) # Create the X,Y matrices for plotting
        Z1 = u[:,:,t] # Pull the Z data matrix for plotting
        Z2 = v[:,:,t] # Pull the Z data matrix for plotting
        fig, (ax1, ax2) = plt.subplots(2,1) # Create the figure with subplots
        # Create the filled countour plot with colormap and manual levels
        cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm,levels=np.arange(plot_u_min,plot_u_max,0.05))
        # cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm)
        fig.colorbar(cf1, ax=ax1) # Add the colorbar
        ax1.set_xlabel('x (mm)') # Label the x axis
        ax1.set_ylabel('y (mm)') # Label the y axis
        utitle = 'Concentration of u at {0:.1f}s'.format(telap)
        ax1.set_title(utitle) # Title the plot
        ax1.set_aspect('equal') # Make the aspect ratio equal
        # Create the filled countour plot with colormap and manual levels
        cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm,levels=np.arange(plot_v_min,plot_v_max,0.05))
        # cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm)
        fig.colorbar(cf2, ax=ax2) # Add the colorbar
        ax2.set_xlabel('x (mm)') # Label the x axis
        ax2.set_ylabel('y (mm)') # Label the y axis
        vtitle = 'Concentration of v at {0:.1f}s'.format(telap)
        ax2.set_title(vtitle) # Title the plot
        ax2.set_aspect('equal') # Make the aspect ratio equal
        plt.subplots_adjust(hspace=0.75,left=-0.5)
    else: # If the Boolean is set to plot only v
        x = np.arange(0,Lx,dx) # Create the x data for plotting
        y = np.arange(0,Ly,dy) # Create the y data for plotting
        X, Y = np.meshgrid(y,x) # Create the X,Y matrices for plotting
        Z = v[:,:,t] # Pull the Z data matrix for plotting
        fig, ax = plt.subplots() # Create the figure with subplots
        # Create the filled countour plot with colormap and manual levels
        cf = ax.contourf(X,Y,Z,cmap=cm.coolwarm,levels=np.arange(plot_v_min,plot_v_max,0.05))
        fig.colorbar(cf, ax=ax) # Add the colorbar
        ax.set_xlabel('x (mm)') # Label the x axis
        ax.set_ylabel('y (mm)') # Label the y axis
        vtitle = 'Concentration of v at {0:.1f}s'.format(telap)
        ax.set_title(vtitle) # Title the plot
        ax.set_aspect('equal') # Make the aspect ratio equal 
    # plt.show() # This shows the plots as the code is running
    fname = 'BZ-TF/Patterns_%06d.png' % t # Create the file name for each plot
    print('Saving frame', fname) # Print the status update
    fig.savefig(fname, dpi=300) # Save the image
    files.append(fname) # Update the filename
    plt.close(fig) # Close the image so it doesn't show while the code is running
    t=t+math.floor(dtWindow/dt)            # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed   

#%% Export images to videos

convert.export("BZ-TF")