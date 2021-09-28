"""
Parameters:
u - Salt
T - Temperature
w - Ice
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
    names = list(range(0, ny))
    for i, n in enumerate(data):
        variableName = dataName[i]
        folder = f'data/{variableName}/data_{t}.csv'          # Specify where to save the array
        data[i] = pd.read_csv(folder, names = names)

#%% Model Parameters

# The time and space parameters for the model go here.
res = 0.5   #0.5    # This sets the resolution for the simulation in mm/step.
Lx = 50     #30     # This sets the x length of the simulation in mm.
Ly = 81     #85     # This sets the y length of the simulation in mm.
tEnd = 75           # This sets the duration of the simulation in s.
dt = 0.05   #0.05   # The time step for the calculation in s.
dtWindow = 0.5    # This sets how often to update the plot in s.

# These are internal parameters for the simulation.
nx = math.floor(Lx/res)             # The x dimension for storage.
ny = math.floor(Ly/res)             # The y dimension for storage.
nt = math.floor(tEnd/dt)+1         # The t dimension for storage.

u = np.zeros((nx,ny))              # Define the array for u, T, w data.
T = np.zeros((nx,ny))              # Define the array for u, T, w data.
w = np.zeros((nx,ny))              # Define the array for u, T, w data.
l = np.zeros((2,nt))               # Define the array for l data.

data = [u,T,w]
dataName = ['u','T','w']

#%% Chemical Parameters

# https://dtrx.de/od/diff/
# Diffusion rate of Water 1.05 to 0.187 μm^2/ms for water in temp of 0.35 °C to –30.65 °C
# However, we assume constant rate of diffusion. Is this accpetable?

diff = [0.002, 0.001, 0] # This sets the diffusion in mm^2/s.
vel = [0, 0, 0]          # This adjusts the advection velocity. 1 for On, 0 for Off

xvel = 0.0    # This set the x advection velocity in mm/s.
yvel = 0.833  # 1.5/18 cm/s = 0.833 mm/s from Harvard Experiment

#%% Initial Conditions
# https://nsidc.org/cryosphere/seaice/characteristics/brine_salinity.html
# The average salinity of the ocean typically varies from 32 to 37 psu
# For every 5 psu increase in salinity, the freezing point decreases by 0.28 degrees Celsius  
        
oceanSalinity = 32          # Initial Salt Concentration in psu
brineSalinity = 60          # Salt concentration of brine is considered 50+ psu

oceanTemperature = 0        # 10 Gallons of water between -2 and 1C
brineTemperature = -20      # Freezer was maintained between −17◦C and −20◦C

ratio = -5/0.28            # Salinity(psu)/Temperature(◦C) ratio

u[:,:] = oceanSalinity      # Set the initial Salt concentration          
T[:,:] = oceanTemperature   # Set the initial Temperature 

#%% Run the Simulation

t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
while telap <= tEnd:

    # Save the Data Arrays
    saveData(t)

    # Get the length data
    l[0,t] = t*dt
    l[1,t] = functions.length(l, data[2],ny)

    # Get the Data for t+1 by calling the pde() function in the file called "functions.py"
    data = functions.pde(data,diff,vel,xvel,yvel,nx,ny,dt,res,telap,brineSalinity,brineTemperature,ratio)
    
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
        
        functions.plot(folder,tEnd, Lx, Ly, res, data[0],data[1],data[2], l, dt, dtWindow, telap)
        t=t+math.floor(dtWindow/dt)            # Increment the storage counter
        telap=t*dt                             # Increment the simulation time elapsed
else:    
    time = int(t/dt)        # Convert time in seconds to index value
    loadData(time)          # Load the Data Arrays
    
    functions.plot(folder,tEnd, Lx, Ly, res, data[0],data[1],data[2], l, dt, dtWindow, t)
    
print("Run Time:", timer()-start)
    
#%% Make Video From Images
#functions.export('images')
    