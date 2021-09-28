"""
So Heres the General Idea:
    You have a big ass tensor of size (100,162,15001) that takes up a lot of ram. 
    If we save a (100,162) every time step, it's is very slow. We need a balance.
Steps:
    1.  Save a 3D array to txt file as a 2D array every n time steps.
        This will be a nth the number of files but more ram than before.
    2.  Write the nth time step to a place holder. Erase the (100,162,n-1) array
        Write the next calculation to the (100,162,0) space
"""

#%% Imports

import math
import numpy as np
import functions
#import pandas as pd

from timeit import default_timer as timer
#start = timer()

#%% Data Handling

def saveData(t):
    start = timer()
    # reshaping the array from 3D matrice to 2D matrice.
    # (5, 4, 3) matirc becomes (5,12)    
    for i, n in enumerate(data):
        arr_reshaped = n.reshape(n.shape[0], -1)
        variableName = dataName[i]
        folder = f'data/{variableName}/data_{t}.csv'          # Specify where to save the array
        np.savetxt(folder, arr_reshaped) # Saves the array to text file
    print("Save Time:", timer()-start)

    
def loadData(t):
    for i, n in enumerate(data):
        variableName = dataName[i]
        folder = f'data/{variableName}/data_{t}.csv'          # Specify where to save the array
        loaded_arr = np.loadtxt(folder)
        data[i] = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // n.shape[2], n.shape[2])
        
def splitNum(nt):
    n = 8
    while nt%n != 0:
        n+=1
    return n
    
#%% Model Parameters
# The time and space parameters for the model go here.
res = 0.5   #0.5    # This sets the resolution for the simulation in mm/step.
Lx = 50     #30     # This sets the x length of the simulation in mm.
Ly = 81     #85     # This sets the y length of the simulation in mm.
tEnd = 75           # This sets the duration of the simulation in s.
dt = 0.005   #0.05   # The time step for the calculation in s.
dtWindow = 0.5    # This sets how often to update the plot in s.

# These are internal parameters for the simulation.
nx = math.floor(Lx/res)             # The x dimension for storage.
ny = math.floor(Ly/res)             # The y dimension for storage.
nt = math.floor(tEnd/dt)+1         # The t dimension for storage.

l = np.zeros((2,nt))                     # Define the array for l data.

split = 100#splitNum(nt)                                # Defines the time size of data arrays
u = np.zeros((nx,ny,split))              # Define the array for u data.
T = np.zeros((nx,ny,split))              # Define the array for T data.
w = np.zeros((nx,ny,split))              # Define the array for w data.

data = [u,T,w]
dataName = ['u','T','w']

#%% Chemical Parameters

diff = [0.002, 0.001, 0] # This sets the diffusion in mm^2/s.
vel = [1, 1, 0]          # This adjusts the advection velocity. 1 for On, 0 for Off

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

ratio = -5/0.28             # Salinity(psu)/Temperature(◦C) ratio

u[:,:,0] = oceanSalinity      # Set the initial Salt concentration          
T[:,:,0] = oceanTemperature   # Set the initial Temperature
#%% Run the Simulation

t = 0        # This sets the first time point for the calculation.
telap = 0    # This sets the time elapsed for the simulation.
n = 0        # This sets what index of the array we are using
fileNum = 0  # This sets which file number to assign to the save data
while telap < tEnd:

    # Get the length data
    l[0,t] = t*dt
    l[1,t] = functions.length(l, data[2][:,:,n],ny)

    # Get the Data for t+1 by calling the pde() function in the file called "functions.py"
    try:
        data[0][:,:,n+1],data[1][:,:,n+1],data[2][:,:,n+1] = functions.pde(data,diff,vel,xvel,yvel,nx,ny,dt,res,telap,brineSalinity,brineTemperature,ratio,n)
    except:
        print("Saving File",fileNum)
        saveData(fileNum)
        fileNum += 1
        # Take the last index and store it in the first
        data[0][:,:,0],data[1][:,:,0],data[2][:,:,0] = functions.pde(data,diff,vel,xvel,yvel,nx,ny,dt,res,telap,brineSalinity,brineTemperature,ratio,n)
        n = -1
        
    
    n+=1              # Increment the index counter
    t+=1              # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
    
    # This displays the step and time being plotted. 
    print('Step: {0} Time: {1:0.3f}s'.format(t,telap))
#%% Plot The Data

#loadData(int(nt/split)-2)
t = tEnd
folder = 'images' # Specifies which folder you would like to save the plots in
data2 = [data[n][:,:,n],data[1][:,:,n],data[2][:,:,n]]
functions.plot(folder,tEnd, Lx, Ly, res, data2, l, dt, dtWindow, t)
#print("Run Time:", timer()-start)
    
