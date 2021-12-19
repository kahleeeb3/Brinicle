"""
Thoughts:
    (1) Diffusion by temp no longer has diffusion barrier due to ice.
        Need to find a way to do both maybe?
    (2) How do I do diffusion by Salt Concentration?
        fyvel now needs to be a matrix: fyvel_matrix = fyvel * up
    (3) Brine Temp Too Cold?
"""
    

#%% Imports

import math
import numpy as np
import PDES as f
import plot as p

#%% Analyzing run-time
#from timeit import default_timer as timer
#start = timer()

#%% Model Parameters

# The time and space parameters for the model go here.
tEnd = 75         # This sets the duration of the simulation in s.
dt = 0.05         # The time step for the calculation in s.
dtWindow = 0.5    # This sets how often to update the plot in s.

param = {
"res": 0.5,         # This sets the resolution for the simulation in mm/step.
"Lx" : 50,          # This sets the x length of the simulation in mm.
"Ly" : 81          # This sets the y length of the simulation in mm.
}

# These are internal parameters for the simulation.
param.update( {
"nx": math.floor(param["Lx"]/param["res"]),  # The x dimension for storage.
"ny": math.floor(param["Ly"]/param["res"]),   # The y dimension for storage.
"nt": int(tEnd/dtWindow),                     # The t dimension for storage.
"ntl": math.floor(tEnd/dt)  # The t dimension for L.
})

S0 = np.zeros((param["nx"],param["ny"]))    # Define the array for S data.
I0 = np.zeros((param["nx"],param["ny"]))    # Define the array for I data.
T0 = np.zeros((param["nx"],param["ny"]))    # Define the array for T data.

L = np.zeros((2,param["ntl"]))                # Define the array for length data.


S = np.zeros((param["nx"],param["ny"],param["nt"]))     # Store Data for S over time.
I = np.zeros((param["nx"],param["ny"],param["nt"]))     # Store Data for I over time.
T = np.zeros((param["nx"],param["ny"],param["nt"]))     # Store Data for T over time.

# Store all Info in a dictonary
data = {
"S":{
     "Name":"Salt",
     "Curr":S0
     },
"I":{
     "Name":"Ice",
     "Curr":I0
     },
"T":{
     "Name":"Temperature",
     "Curr":T0
     }
}      # Dictionary of variables

del S0,I0,T0 # remove values from memory

#%% Initial Conditions
# https://nsidc.org/cryosphere/seaice/characteristics/brine_salinity.html
# The average salinity of the ocean typically varies from 32 to 37 psu
# For every 5 psu increase in salinity, the freezing point decreases by 0.28 degrees Celsius  
param.update({    
"oceanSalinity": 32,          # Initial Salt Concentration in psu
"brineSalinity": 60,          # Salt concentration of brine is considered 50+ psu
"oceanTemperature": 0,        # 10 Gallons of water between -2 and 1C
"brineTemperature": -10,     # Freezer was maintained between −17◦C and −20◦C
"ratio": -5/0.28             # Salinity(psu)/Temperature(◦C) ratio
}) 

#%% Chemical Parameters

data["S"]["diff"] = 0.002   # This sets the diffusion in mm^2/s for S
data["I"]["diff"] = 0.000   # This sets the diffusion in mm^2/s for I
data["T"]["diff"] = 0.001   # This sets the diffusion in mm^2/s for T

# This adjusts the advection velocity. 1 for On, 0 for Off
data["S"]["vel"] = 1
data["I"]["vel"] = 0
data["T"]["vel"] = 1        

param["xvel"] = 0.0  # This set the x advection velocity in mm/s.  
param["yvel"] = 0.7    # 1.5/18 cm/s = 0.833 mm/s from Harvard Experiment
#%% Run the Simulation
t = 0       # Index of the length array
telap = 0   # Time elapsed during the simulation.
n = 0       # Frame number

while telap < tEnd:
    # Get the length data
    L[0,t] = telap    # Store time
    L[1,t] = f.length(data["I"]["Curr"],param["ny"]) # Store Len
    
    # if we are at a frame we would like to plot, store the frame
    if telap%dtWindow == 0:
        print(f'Storing Frame at {telap}s')
        S[:,:,n] = data["S"]["Curr"]
        I[:,:,n] = data["I"]["Curr"]
        T[:,:,n] = data["T"]["Curr"]
        n+=1    # increment the frame counter
        
    f._global(data,param,dt) # Calculate the next time step
    
    t+=1                       # Increment the length index
    telap=t*dt                 # Increment the simulation time elapsed

#%% Plot The Data
t = tEnd    # This sets the first time point for the Plotting. 
#t = 0       # Set to 0 to produce all images

folder = 'images' # Specifies which folder you would like to save the plots in
if t == tEnd:
    data2 = [data["S"]["Curr"],data["T"]["Curr"],data["I"]["Curr"],L]
    time = [ tEnd , dt , dtWindow, telap]
    p.plot(folder,data2,param,time)
else:
    telap = 0    # This sets the time elapsed for the simulation.
    n = 0
    while telap < tEnd:
        data2 = [S[:,:,n],T[:,:,n],I[:,:,n],L]
        
        n+=1                # Increment the time step        
        t+=1                # Increment the storage counter
        telap=t*dtWindow    # Increment the simulation time elapsed
        
        time = [ tEnd , dt , dtWindow, telap]
        p.plot(folder,data2,param,time)
        

# Make Video
p.export('images')
#Evaluate Run-Time
#print("Run Time:", timer()-start)
