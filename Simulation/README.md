# Python Simulations
Note: * I found the [Syder IDE](https://www.spyder-ide.org/) to work well for this project as it allows the user access to a variable explorer to better analyze the values in our arrays.*
## 1. Mathematical Model
Our simulation works by implementing Newton's Method of iteration. We use an initial condition of our function to estimate the value of our function over a change in time such that: 
<!-----------LATEX IN HTML----------->
<div align ="center"> 
    <img src="https://latex.codecogs.com/gif.latex?\dpi{150}&space;\bg_black&space;\fn_jvn&space;\boxed{x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}.}"/>
</div>
<!--------------------------------------->

## 2. Model Parameters
The 3 parameter of the function are defined as `S,I,& T` which represent the concentration of `Salt, Ice, and Temperature` respectively. The value of these parameters at a given time step are stored as second rank tensors whose size is defined by the resolution of our image, and the length and width of our system:
```python
S0 = np.zeros((nx,ny))              # Define the array for S data.
I0 = np.zeros((nx,ny))              # Define the array for I data.
T0 = np.zeros((nx,ny))              # Define the array for T data.
```
After each calculation, these parameters are stored in their respective index of a third rand whose third dimension size is defined by our simulation time in seconds divided by how often we would like to plot.
```python
S = np.zeros((nx,ny,nt))         # Store Data for S over time.
I = np.zeros((nx,ny,nt))         # Store Data for I over time.
T = np.zeros((nx,ny,nt))         # Store Data for T over time.
```
Note on RAM usage: By using `S0,I0,& T0` solely for calculation at each time step, and storing only the frames we would like to plot in `S,I,& T`, we enables the simulation to run at higher resolutions due to the massive cut in RAM usage incurred by defining three individual third rank tensors which store all time step calculations in a given parameter.
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