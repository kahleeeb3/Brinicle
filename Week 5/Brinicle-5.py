"""
Parameters:
u - Salt
v - Brine
w - Ice
T - Temperature
"""

#%% Imports

import math
import numpy as np
import functions
import pandas as pd

from timeit import default_timer as timer

start = timer()

#%% Data Handling

def saveData(t):
    for i, n in enumerate(data):
        variableName = dataName[i]
        folder = f'data/{variableName}/data_{t}.csv'          # Specify where to save the array
        np.savetxt(folder, n, delimiter=',')                    # Save the array to a csv.

def loadData(t):
    for i, n in enumerate(data):
        variableName = dataName[i]
        folder = f'data/{variableName}/data_{t}.csv'          # Specify where to save the array
        data[i] = pd.read_csv(folder)

#%% Model Parameters

# The time and space parameters for the model go here.
res = 0.5   #0.5    # This sets the resolution for the simulation in mm/step.
Lx = 50     #30     # This sets the x length of the simulation in mm.
Ly = 81     #85     # This sets the y length of the simulation in mm.
tEnd = 30           # This sets the duration of the simulation in s.
dt = 0.05   #0.05   # The time step for the calculation in s.
dtWindow = 0.5    # This sets how often to update the plot in s.

# These are internal parameters for the simulation.
dx = res                           # The x resolution for the calculation.
dy = res                           # The y resolution for the calculation.
nx = math.floor(Lx/dx)             # The x dimension for storage.
ny = math.floor(Ly/dy)             # The y dimension for storage.
nt = math.floor(tEnd/dt)+1         # The t dimension for storage.

l = np.zeros((2,nt))                 # Define the array for l data.

data = [np.zeros((nx+1,ny)),np.zeros((nx+1,ny)),np.zeros((nx+1,ny))]
dataName = ['u','v','w']

#%% Chemical Parameters

diff = [0.002,0.002,0] # This sets the diffusion in mm^2/s.
vel = [0.0,0.7,0.0]   # This adjusts the advection velocity.

xvel = 0.0    # This set the x advection velocity in mm/s.
yvel = 0.5    # This set the y advection velocity in mm/s.

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

data[0][:,:] = 1 #u0 + ampu * np.random.rand(nx,ny) # Adds noise to IC

#%% Run the Simulation

t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap <= tEnd:

    # Save the Data Arrays
    saveData(t)

    # Get the length data
    l[0,t] = t*dt
    l[1,t] = functions.length(l, data[2],res)

    # Get the Data for t+1 by calling the pde() function in the file called "functions.py"
    data = functions.pde(data,diff,vel,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d)
    
    t=t+1              # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    
    # This displays the step and time being plotted. 
    print('Step: {0} Time: {1:0.3f}s'.format(t,telap))

#%% Plot The Data

# Specifies which time plot you would like to make
t = tEnd  # Set this value to 0 to make all plots
folder = 'images' # Specifies which folder you would like to save the plots in

if t == 0:
    telap = 0    # This sets the time elapsed for the simulation.
    while telap <= tEnd:  
        
        loadData(t)          # Load the Data Arrays
        
        functions.plot(folder,tEnd, Lx, Ly, dx, dy, data[0],data[1],data[2], l, dt, dtWindow, telap)
        t=t+math.floor(dtWindow/dt)            # Increment the storage counter
        telap=t*dt                             # Increment the simulation time elapsed
else:    
    time = int(t/dt)        # Convert time in seconds to index value
    loadData(time)          # Load the Data Arrays
    
    functions.plot(folder,tEnd, Lx, Ly, dx, dy, data[0],data[1],data[2], l, dt, dtWindow, t)
    
print("Run Time:", timer()-start)
    
#%% Make Video From Images

#functions.export("images")
    