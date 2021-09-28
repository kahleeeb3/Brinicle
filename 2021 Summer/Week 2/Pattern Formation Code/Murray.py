#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  24 2021

A discretized point model for simulating Turing patterns
Pattern formation model from J. D. Murray
https://www.sciencedirect.com/science/article/pii/S0022039602001560
https://core.ac.uk/download/pdf/82316523.pdf

@author: tompkinn

Cat Furr Patterns
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
def PatternsSteadyState(x):
    # Returns u,v for no diffusion.
    up = x[0]
    vp = x[1]
    # These three lines implement the model.
    h=rho*up*vp/(1+up+K*up**2);
    upp = gamma*(a-up-h)
    vpp = alpha*(b-vp)-h
    y = [upp,vpp]
    
    return y

# %% This is the reaction diffusion solver.
def pde(up,vp,gamma,a,b,alpha,rho,K,du,dv,dt,dx,dy):
# This numerically integrates using finite differences and the forward
# Euler method.

    # Set up periodic bc's in y.
    # These pad the array with the last items from the other side,
    # Thus the last item wraps around to the other side.
    upy = np.column_stack((up[:,-1],up,up[:,1]))
    vpy = np.column_stack((vp[:,-1],vp,vp[:,1]))
    # Make boundaries periodic in x.
#    upxy = np.vstack((upy[-1,:], upy, upy[1,:]))
#    vpxy = np.vstack((vpy[-1,:], vpy, vpy[1,:]))
        
    # Set up zero flux bc's in y.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
#    upy = np.column_stack((up[:,2],up,up[:,-2]))
#    vpy = np.column_stack((vp[:,2],vp,vp[:,-2]))
    # Make boundaries zero flux in x.
    upxy = np.vstack((upy[2,:], upy, upy[-2,:]))
    vpxy = np.vstack((vpy[2,:], vpy, vpy[-2,:]))
            
    # Perform finite differences. (Diffusion)
    # On axis terms for u
    # Calculate differences
    uxx=(1/dx**2)*(upxy[:,2:]+upxy[:,0:-2]-2*upxy[:,1:-1])
    uxx=uxx[1:-1,:] # Remove extra rows
    # Calculate differences
    uyy=(1/dy**2)*(upxy[2:,:]+upxy[0:-2,:]-2*upxy[1:-1,:])
    uyy=uyy[:,1:-1] # Remove extra columns
    # Diagonal terms for u
    # Calculate differences
    uxy=(1/(dx**2+dy**2))*(upxy[2:,2:]+upxy[0:-2,0:-2]-2*upxy[1:-1,1:-1])
    uyx=(1/(dx**2+dy**2))*(upxy[0:-2,2:]+upxy[2:,0:-2]-2*upxy[1:-1,1:-1])
        
    # On axis terms for v
    # Calculate differences
    vxx=(1/dx**2)*(vpxy[:,2:]+vpxy[:,0:-2]-2*vpxy[:,1:-1])
    vxx=vxx[1:-1,:] # Remove extra rows
    # Calculate differences
    vyy=(1/dy**2)*(vpxy[2:,:]+vpxy[0:-2,:]-2*vpxy[1:-1,:])
    vyy=vyy[:,1:-1] # Remove extra columns
    # Diagonal terms for v
    # Calculate differences
    vxy=(1/(dx**2+dy**2))*(vpxy[2:,2:]+vpxy[0:-2,0:-2]-2*vpxy[1:-1,1:-1])
    vyx=(1/(dx**2+dy**2))*(vpxy[0:-2,2:]+vpxy[2:,0:-2]-2*vpxy[1:-1,1:-1])
            
    # Combine diffusion with kinetics using the dangerous forward
    # Euler algorithm.
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15;  # Set ff=1 to turn the fudge-factor off.
    h=rho*up*vp/(1+up+K*up**2);
    up=up+dt*(du*(uxx+uyy+(uxy+uyx)/ff)) + dt*(gamma*(a-up-h))
    vp=vp+dt*(dv*(vxx+vyy+(vxy+vyx)/ff)) + dt*(alpha*(b-vp)-h)
        
    return [up, vp]

# %% Model Parameters

# The time and space parameters for the model go here.
res = 0.05     # This sets the resolution for the simulation in mm/step.
Lx = 5         # This sets the x length of the simulation in mm.
Ly = 12         # This sets the y length of the simulation in mm.
tEnd = 10     # This sets the duration of the simulation in s.
dt = 0.005       # The time step for the calculation in s.
dtWindow = 0.1 # This sets how often to update the plot in s.

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
Relevant Parameters to change
"""
a=92        # 92
b=64        # 64
d=15        # 15    # ratio between diffusion rates

# The chemistry parameters of the model go here.
"""
Changes the spacing
"""
du = 0.002    # This sets the diffusion for u in mm^2/s.
dv = d*du     # This sets the diffusion for v in mm^2/s.

#xvel = 0.0      # This set the x advection velocity in mm/s.
#yvel = 0.0      # This set the y advection velocity in mm/s.

# Kinetic parameters # Defaults
alpha=1.5   # 1.5
K=0.1       # 0.1
rho=18.5    # 18.5
gamma=9     # 9

# %% Steady State Solution

# The code starts by finding the steady state values in the absence of
# diffusion. This is the fixed point for the system.
initGuess = [1,1]
uv0 = fsolve(PatternsSteadyState,initGuess)
u0 = uv0[0]   # Steady state value for u
v0 = uv0[1]   # Steady state value for v

# %% Initial Conditions

# This sets the initial conditions for the simulation. Initially every spot
# is set to the same values with some optional random variational
# of amplitude amp.
# u0 = 0.3        # Arbitrary initial value
# v0 = 0.3        # Arbitrary initial value
"""
Heterogeneity is on
If you make everything perfectly homogenous, you will not get patterns.
"""
u_amp = u0/3      # Set to zero to start everything the same
v_amp = v0/3      # Set to zero to start everything the same
u[:,:,0] = u0 + u_amp * np.random.rand(nx,ny)
v[:,:,0] = v0 + v_amp * np.random.rand(nx,ny)
uC = u0/4          # The u value for the different spots.
vC = v0/4          # The v value for the different spots.

"""
Change the shape and size
"""

#u[:,:,0] = convert.image("image.png",nx,ny,180)/255 * u0
#v[:,:,0] = convert.image("image.png",nx,ny,90)

# This is the middle spot of specified width
spotwidth=10 # This is half width
width = 90
spotleft=int(np.round(nx/2))-spotwidth - width # Determine the left edge
spotright=int(np.round(nx/2))+spotwidth + width # Determine the right edge
spottop=int(np.round(ny/2))-spotwidth # Determine the top edge
spotbottom=int(np.round(ny/2))+spotwidth # Determine the bottom edge

#u[spotleft:spotright,spottop:spotbottom,0]=uC # Create the initial spot in u
#v[spotleft:spotright,spottop:spotbottom,0]=vC # Create the initial spot in y

# %% This section runs the simulation itself.
t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap < tEnd:    
       
    # From here on down is the nuts and bolts of the simulation.
    # Update u and v for the next time output point defined by nWindow
    # using the pde function defined at the beginning of this file
    [u[:,:,t+1],v[:,:,t+1]]=pde(u[:,:,t],v[:,:,t],gamma,a,b,alpha,rho,K,du,dv,dt,dx,dy)
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
plot_u_min = -10
plot_u_max = 10
plot_v_min = -10
plot_v_max = 10

text = f'a = {a}\nb = {b}\nd = {d}'
textx,texty,size = 17,0,11


while telap <= tEnd:

    # Make the plots
    if displayUandV: # If the Boolean is set to plot both u and v
        x = np.arange(0,Lx,dx) # Create the x data for plotting
        y = np.arange(0,Ly,dy) # Create the y data for plotting
        X, Y = np.meshgrid(y,x) # Create the X,Y matrices for plotting
        Z1 = u[:,:,t] # Pull the Z data matrix for plotting
        Z2 = v[:,:,t] # Pull the Z data matrix for plotting
        fig, (ax1, ax2) = plt.subplots(2,1) # Create the figure with subplots
        # Create the filled countour plot with colormap and manual levels
        # cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm,levels=np.arange(plot_u_min,plot_u_max,0.05))
        cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm)
        fig.colorbar(cf1, ax=ax1) # Add the colorbar
        ax1.set_xlabel('x (mm)') # Label the x axis
        ax1.set_ylabel('y (mm)') # Label the y axis
        utitle = 'Concentration of u at {0:.1f}s'.format(telap)
        ax1.set_title(utitle) # Title the plot
        ax1.set_aspect('equal') # Make the aspect ratio equal
        # Create the filled countour plot with colormap and manual levels
        # cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm,levels=np.arange(plot_v_min,plot_v_max,0.05))
        cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm)
        fig.colorbar(cf2, ax=ax2) # Add the colorbar
        ax2.set_xlabel('x (mm)') # Label the x axis
        ax2.set_ylabel('y (mm)') # Label the y axis
        vtitle = 'Concentration of v at {0:.1f}s'.format(telap)
        ax2.set_title(vtitle) # Title the plot
        ax2.set_aspect('equal') # Make the aspect ratio equal
        plt.subplots_adjust(hspace=0.75,left=-0.5)
        plt.text(textx, texty, text, fontsize=size)
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
    fname = 'Murray/Patterns_%06d.png' % t # Create the file name for each plot
    print('Saving frame', fname) # Print the status update
    fig.savefig(fname, dpi=300) # Save the image
    files.append(fname) # Update the filename
    plt.close(fig) # Close the image so it doesn't show while the code is running
    t=t+math.floor(dtWindow/dt)            # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    
#%% Export images to videos

convert.export("Murray")
