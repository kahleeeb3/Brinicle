#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parameters:
u - Salt
v - Brine
w - Ice

The equations then become:
du/dt = a*u*v 
dv/dt = -b*u*v+c*u*v
dw/dt = d*u*v
"""

# %% Imports

import math
import numpy as np
#from scipy.optimize import fsolve
#import matplotlib.pyplot as plt
#from matplotlib import cm
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
    vyvel=dt*(vel_v*yvel/dy)*(vpxy[:,2:]-vpxy[:,1:-1])
    vyvel=vyvel[1:-1,:] # Remove extra rows
    wyvel=dt*(vel_w*yvel/dy)*(wpxy[:,0:-2]-wpxy[:,1:-1])
    wyvel=wyvel[1:-1,:] # Remove extra rows
            
    # Combine diffusion with advection and kinetics using
    # the forward Euler algorithm.
    """
    Change:
        du/dt = -k1*u*v 
        dv/dt = -k1*u*v
        dw/dt = +k2*u*v
    
    The equations then become:
        du/dt = -a*u*v 
        dv/dt = -b*u*v+c*u*v
        dw/dt = d*u*v
    """
    up=up+(uxx+uyy+uxy+uyx)+(uxvel+uyvel) + dt*(-a*up*vp)     #(-k1*up*vp)
    vp=vp+(vxx+vyy+vxy+vyx)+(vxvel+vyvel) + dt*(-b*up*vp+c*up*vp)     #(-k1*up*vp)
    wp=wp+(wxx+wyy+wxy+wyx)+(wxvel+wyvel) + dt*(d*up*vp)     #(+k1*up*vp)
    
    # Here you apply any constant conditions
    # Holding the top left edge at a constant value
    #tl_u = 0.0
    tl_v = 0.0
    # From fraction to end (top)
    """
    up[0,int(np.round(5*ny/10)):] = tl_u
    """
    spotwidth=2 # This is half width in steps
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    vp[spotleft:spotright,ny-5:ny-2]= 1
    
    up[0,int(np.round(5*ny/10)):] = tl_v
    # Holding the bottom left edge at a constant value
    bl_u = 0.0
    bl_v = 0.0
    # From begining (bottom) to fraction
    up[0,0:int(np.round(5*ny/10))] = bl_u
    vp[0,0:int(np.round(5*ny/10))] = bl_v
        
    return [up, vp, wp]
    

# %% Model Parameters

# The time and space parameters for the model go here.
res = 0.5        # This sets the resolution for the simulation in mm/step.
Lx = 50         # This sets the x length of the simulation in mm.
Ly = 81          # This sets the y length of the simulation in mm.
tEnd = 75   #123      # This sets the duration of the simulation in s.
dt = 0.005         # The time step for the calculation in s.
dtWindow = 0.5    # This sets how often to update the plot in s.

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
"""
Parameters:
u - Salt
v - Brine
w - Ice
"""
du = 0.002    # This sets the diffusion for u in mm^2/s.
dv = 0.002    # This sets the diffusion for v in mm^2/s.
dw = 0        # This sets the diffusion for v in mm^2/s.
xvel = 0.0    # This set the x advection velocity in mm/s.
yvel = 0.02    # This set the y advection velocity in mm/s.
vel_u = 0     # This adjusts the advection velocity for u.
vel_v = 1     # This adjusts the advection velocity for v.
vel_w = 0     # This adjusts the advection velocity for w.

"""
Changing the rates:
    du/dt = a*u*v 
    dv/dt = -b*u*v+c*u*v
    dw/dt = d*u*v
"""
# 
a = 1.0 # ?.?   # consumption of salt water
b = 0.0 # 0.0   # loss of brine
c = 1.0 # ?.?   # production of brine
d = 0.8 # 0.8   # production of ice per liter


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
"""
Parameters:
u - Salt
v - Brine
w - Ice
"""

u[:,:,0] = 0 #u0 + ampu * np.random.rand(nx,ny) # Adds noise to IC
v[:,:,0] = 0 #v0 + ampv * np.random.rand(nx,ny) # Adds noise to IC
w[:,:,0] = 0 #w0 + ampw * np.random.rand(nx,ny) # Adds noise to IC
uC = 0.5      # The u value for the different spots.
vC = 0.5      # The v value for the different spots.
wC = 0.5      # The v value for the different spots.

# This is the middle spot of specified width
w[0:nx,ny-2:ny,0]=1 # Create the initial spot in u

spotwidth=2 # This is half width in steps
spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
#v[spotleft:spotright,ny-5:ny-2,0]=2 # Create the initial spot in v
u[:,:,0] = 1

#spottop=int(np.round(ny/2))-spotwidth    # Determine the top edge
#spotbottom=int(np.round(ny/2))+spotwidth # Determine the bottom edge

# data for the length
l = np.zeros((2,nt-1))


# %% This section runs the simulation itself.
t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap < tEnd:    
       
    # From here on down is the nuts and bolts of the simulation.
    # Update u, v, and w for the next time output
    # using the pde function defined at the beginning of this file
    [u[:,:,t+1],v[:,:,t+1],w[:,:,t+1]]=ard_pde(u[:,:,t],v[:,:,t],w[:,:,t],du,dv,dw,xvel,yvel,dt,dx,dy)
    # calc the length of the thing using the length function
    l[0,t] = t*dt
    l[1,t] = convert.improvedlength(l, w[:,:,t])
    ups=u[1:4,1:4,t+1] # Samples for debugging
    vps=v[1:4,1:4,t+1] # Samples for debugging
    wps=w[1:4,1:4,t+1] # Samples for debugging
    t=t+1              # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    # This displays the step and time being plotted. 
    print('Step: {0} Time: {1:0.3f}s'.format(t,telap))
        

#%% Export images to videos

convert.plot("images",tEnd, np,Lx, Ly, dx, dy, u, v , w, l, dt, dtWindow)

#%% Export images to videos

convert.export("images")
