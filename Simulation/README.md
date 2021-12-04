# Python Simulations
The general functions of the code is to define 3 parameter (u,v,w) for the Salt, Temperature, and Ice Concentration. These are stored as 3D arrays that are defined in the code below
## Model Parameters
```python
# The time and space parameters for the model go here.
res = 0.5 # This sets the resolution for the simulation in mm/step.
Lx = 50 # This sets the x length of the simulation in mm.
Ly = 81 # This sets the y length of the simulation in mm.
tEnd = 75 # This sets the duration of the simulation in s.
dt = 0.05 # The time step for the calculation in s.
dtWindow = 0.5 # This sets how often to update the plot in s.

# These are internal parameters for the simulation.
nx = math.floor(Lx/res) # The x dimension for array.
ny = math.floor(Ly/res) # The y dimension for array.
nt = math.floor(tEnd/dt) # The t dimension for array.

l = np.zeros((2,nt)) # Define the array for length data.

split = int(tEnd/dtWindow) # Defines the time size of data arrays
u = np.zeros((nx,ny,split)) # Define the array for u data.
T = np.zeros((nx,ny,split)) # Define the array for T data.
w = np.zeros((nx,ny,split)) # Define the array for w data.

u0 = np.zeros((nx,ny)) # Define the initial array for u data.
T0 = np.zeros((nx,ny)) # Define the initial array for u data.
w0 = np.zeros((nx,ny)) # Define the initial array for u data.
```
## Chemical Parameters
```python
diff = [0.002, 0.001, 0] # This sets the diffusion in mm^2/s.
vel = [1, 1, 0]          # This adjusts the advection velocity. 1 for On, 0 for Off

xvel = 0.0    # This set the x advection velocity in mm/s.
yvel = 0.7    # 1.5/18 cm/s = 0.833 mm/s from Harvard Experiment
```
## Initial Conditions
```python      
oceanSalinity = 32 # Initial Salt Concentration in psu
brineSalinity = 60 # Salt concentration of brine is considered 50+ psu

oceanTemperature = 0 # 10 Gallons of water between -2 and 1C
brineTemperature = -20 # Freezer was maintained between −17◦C and −20◦C

ratio = -5/0.28 # Salinity(psu)/Temperature(◦C) ratio
```