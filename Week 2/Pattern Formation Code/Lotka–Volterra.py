#%% Lotka–Volterra
"""
The reaction implemented here is the Lotka–Volterra
For convenience the two species x and y are renamed u and v
https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations

The equations then become:
du/du = a * u - b * u * v
dv/dt = d * u * v - g * v

and the fixed point is:
u = g/d
v = a/b

For oscillations use:
a = 2/3
b = 4/3
d = g = 1
"""
# %% Imports

import math
import numpy as np
#from scipy.optimize import fsolve
#import matplotlib.pyplot as plt
#from matplotlib import cm
import convert

# %% This is the reaction diffusion solver.
def pde(up,vp,du,dv,xvel,yvel,dt,dx,dy):
# This numerically integrates using finite differences and the forward
# Euler method.

    # Set up periodic bc's in y.
    # These pad the array with the last items from the other side,
    # Thus the last item wraps around to the other side.
    upy = np.column_stack((up[:,-1],up,up[:,0]))
    vpy = np.column_stack((vp[:,-1],vp,vp[:,0]))
    # Make boundaries periodic in x.
    upxy = np.vstack((upy[-1,:], upy, upy[0,:]))
    vpxy = np.vstack((vpy[-1,:], vpy, vpy[0,:]))
        
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
    """
    The equations then become:
    du/du = a * u - b * u * v
    dv/dt = d * u * v - g * v
    """
    up=up+(uxx+uyy+uxy+uyx)+(uxvel+uyvel) + dt*(a * up - b * up * vp)#(a + (up**2)*vp - b*up - up)
    vp=vp+(vxx+vyy+vxy+vyx)+(vxvel+vyvel) + dt*(d * up * vp - g * vp)#(b*up - (up**2)*vp)
    
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
res = 0.05        # This sets the resolution for the simulation in mm/step.
Lx = 5         # This sets the x length of the simulation in mm.
Ly = 12         # This sets the y length of the simulation in mm.
tEnd = 7       # This sets the duration of the simulation in s.
dt = 0.05        # The time step for the calculation in s.
dtWindow = 0.05   # This sets how often to update the plot in s.

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
du = 1*0.002    # This sets the diffusion for u in mm^2/s.
dv = 1*du       # This sets the diffusion for v in mm^2/s.
xvel = 0.0      # This set the x advection velocity in mm/s.
yvel = 0.0      # This set the y advection velocity in mm/s.
a = 2/3
b = 4/3
d = g = 1

# %% Steady State Solution

# The code starts by finding the steady state values in the absence of
# diffusion. This is the fixed point for the system.
initGuess = [1,1]
#uv0 = fsolve(BrusselatorSteadyState,initGuess)
"""
and the fixed point is:
u = g/d
v = a/b
"""
u0 = g/d   #uv0[0]   # Steady state value for u
v0 = a/b   #uv0[1]   # Steady state value for v

# %% Initial Conditions

# This sets the initial conditions for the simulation. Initially every spot
# is set to the same values with some optional random variational
# of amplitude amp.
# u0 = 0.3        # Arbitrary initial value
# v0 = 0.3        # Arbitrary initial value
"""
Anything Non-Zero (any heterogeneity at all) will have no patterns
Must be perfectly homogenous
"""
u_amp =u0/3*0      # Set to zero to start everything the same
v_amp =v0/3*0      # Set to zero to start everything the same
u[:,:,0] = u0 + u_amp * np.random.rand(nx,ny)
v[:,:,0] = v0*0 + v_amp * np.random.rand(nx,ny)
uC = u0/4          # The u value for the different spots.
vC = v0/4          # The v value for the different spots.

#u[:,:,0] = 2
#v[:,:,0] = 3
# This is the middle spot of specified width
spotwidth=20 # This is half width in steps
spotleft=int(np.round(nx/2))-spotwidth # Determine the left edge
spotright=int(np.round(nx/2))+spotwidth # Determine the right edge
spottop=int(np.round(ny/2))-spotwidth # Determine the top edge
spotbottom=int(np.round(ny/2))+spotwidth # Determine the bottom edge
#u[spotleft:spotright,spottop:spotbottom,0]=uC # Create the initial spot in u
v[spotleft:spotright,spottop:spotbottom,0]=1 # Create the initial spot in y

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
    
#%% Export images to videos

convert.plot("LotkaVolterra",tEnd, displayUandV, np, Lx, Ly, dx, dy, u, v, dt, dtWindow)

#%% Export images to videos

convert.export("LotkaVolterra")
