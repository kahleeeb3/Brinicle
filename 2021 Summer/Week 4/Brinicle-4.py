"""
Parameters:
u - Salt
v - Brine
w - Ice

The partial differential equations:
du/dt = a*u*v 
dv/dt = -b*u*v+c*u*v
dw/dt = d*u*v
"""

#%% Imports

import math
import numpy as np
import functions
#from numba import jit, cuda

#%% Model Parameters

# The time and space parameters for the model go here.
res = 0.5   #0.5    # This sets the resolution for the simulation in mm/step.
Lx = 50     #30     # This sets the x length of the simulation in mm.
Ly = 81     #85     # This sets the y length of the simulation in mm.
tEnd = 75           # This sets the duration of the simulation in s.
dt = 0.05   #0.05   # The time step for the calculation in s.
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
l = np.zeros((2,nt-1))             # Define the array for l data.

#%% Chemical Parameters

du = 0.002    # This sets the diffusion for u in mm^2/s.
dv = 0.002    # This sets the diffusion for v in mm^2/s.
dw = 0        # This sets the diffusion for v in mm^2/s.

xvel = 0.0    # This set the x advection velocity in mm/s.
yvel = 0.5    # This set the y advection velocity in mm/s.

vel_u = 0.0   # This adjusts the advection velocity for u.
vel_v = 0.7   # This adjusts the advection velocity for v.
vel_w = 0.0   # This adjusts the advection velocity for w.

a = 1.0   # consumption of salt water
b = 0.2   # loss of brine
c = 1.0   # production of brine
d = 1.0   # production of ice per liter


#%% Initial Conditions

# This sets the initial conditions for the simulation. Initially every spot
# is set to the same values with some optional random variational of amplitude.
# This has been commmented out to remove variation.
#u0 = 0        # Arbitrary initial value
#v0 = 0        # Arbitrary initial value
#w0 = 0        # Arbitrary initial value
#ampu = 0      # Set to zero to start everything the same
#ampv = 0      # Set to zero to start everything the same
#ampw = 0      # Set to zero to start everything the same

u[:,:,0] = 1 #u0 + ampu * np.random.rand(nx,ny) # Adds noise to IC
v[:,:,0] = 0 #v0 + ampv * np.random.rand(nx,ny) # Adds noise to IC
w[:,:,0] = 0 #w0 + ampw * np.random.rand(nx,ny) # Adds noise to IC
#uC = 0.5      # The u value for the different spots.
#vC = 0.5      # The v value for the different spots.
#wC = 0.5      # The v value for the different spots.

"""
# Sheet of Ice
w[0:nx,ny-2:ny,0]=1 # Create the initial spot in u

# Brine Source
spotwidth=2 # This is half width in steps
spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
"""

#%% Run the Simulation
t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap < tEnd:    
    
    # Applies the PDE to the array by calling the pde() function in the file called "functions.py"
    [u[:,:,t+1],v[:,:,t+1],w[:,:,t+1]]=functions.pde(u[:,:,t],v[:,:,t],w[:,:,t],du,dv,dw,vel_u,vel_v,vel_w,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d)
    
    # Calculates the length of w using the length function in the file called "functions.py"
    l[0,t] = t*dt
    l[1,t] = functions.length(l, w[:,:,t],res)
    
    ups=u[1:4,1:4,t+1] # Samples for debugging
    vps=v[1:4,1:4,t+1] # Samples for debugging
    wps=w[1:4,1:4,t+1] # Samples for debugging
    
    t=t+1              # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    
    # This displays the step and time being plotted. 
    print('Step: {0} Time: {1:0.3f}s'.format(t,telap))        

#%% Create Plots

# original Plot
# Specifies the time you want plotted. 
plotTime = tEnd    # To make all plots, set this value to 0
functions.plot("images",tEnd, Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow,plotTime)

#functions.plotnew("images", tEnd, np, Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow)

#%% Export images to videos

# Creates a video from images in desired folder
functions.export("images")