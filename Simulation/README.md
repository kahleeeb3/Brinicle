# Python Simulations
Notes:
- I found the [Syder IDE](https://www.spyder-ide.org/) to work well for this project as it allows the user access to a variable explorer to better analyze the values in our arrays.
- The concepts of python dictionaries, global variables, and numpy arrays are essential to understanding the functionality of this simulation
<hr>

## 1. Defining Variables (main.py)
The three variable of the simulation are defined as `S,I,T` which represent the concentration of `Salt, Ice, and Temperature` respectively. The values of these variables at a given time step are stored as second rank tensors `(S0,I0,T0)` whose size is defined by the resolution of our image, and the length and width of our system. We also define a fourth variable, `L`, which stores the length of our brinicle at each time step:
```python
S0 = np.zeros((param["nx"],param["ny"]))    # Define the array for S data.
I0 = np.zeros((param["nx"],param["ny"]))    # Define the array for I data.
T0 = np.zeros((param["nx"],param["ny"]))    # Define the array for T data.
L = np.zeros((2,param["ntl"]))              # Define the array for length data.
```
The values in these variables are stored in a single python dictionary defined as `data`. This dictionary is then passed into a series of functions which give us the value of the variables at a the next time step. This can be simply explained using `Newton's Method`:
<!-----------LATEX IN HTML----------->
<div align ="center"> 
    <img src="https://latex.codecogs.com/gif.latex?\dpi{125}&space;\bg_black&space;\fn_jvn&space;\boxed{x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}.}"/>
</div>
<!-------------------------------------->

When the time step we are calculating is a multiple of `dtWindow` (the plot interval time), these parameters are stored in their respective index of a third rank tensor `(S,I,T)` whose third dimension size `(nt)` is defined by our simulation time in seconds divided by how often we would like to plot.
```python
S = np.zeros((param["nx"],param["ny"],param["nt"]))     # Store Data for S over time.
I = np.zeros((param["nx"],param["ny"],param["nt"]))     # Store Data for I over time.
T = np.zeros((param["nx"],param["ny"],param["nt"]))     # Store Data for T over time.
```
<hr>

**Note on the `param` variable:** *I chose to place all constant values not related to time in a single python dictionary defined as `param`. The simulation can be done without this (as done in the `Summer 2021` folder), however, this cuts down the number of parameters to pass into the `PDES` modules which makes the functions look cleaner and slightly cuts down the run-time.*

**Note on RAM usage:** *By using `S0,I0,& T0` solely for calculation at each time step, and storing only the frames we would like to plot in `S,I,& T`, we enables the simulation to run much faster due to the massive cut in RAM usage incurred by defining three individual third rank tensors which store all time step calculations in a given parameter.*
<hr>

## 2. Running the Simulation
After each calculation, the `S0,I0,T0` values stored in the `data` dictionary are overwritten and used to calculate the next time step:
<!-----------LATEX IN HTML----------->
<div align ="center"> 
    <img src="https://latex.codecogs.com/gif.latex?\dpi{125}&space;\bg_black&space;\fn_jvn&space;\boxed{x_n=x_{n+1}.}"/>
</div>
<!-------------------------------------->

As noted above, this greatly reduces the ram usage of the simulation. After defining the variables and constant values, we can produce a single while-loop at the bottom of our `main.py` file:
```python
#%% Run the Simulation
t = 0       # Index of the length array
telap = 0   # Time elapsed during the simulation.
n = 0       # Frame number
while telap < tEnd:

    # Get the length data
    L[0,t] = telap    # Store time value
    L[1,t] = f.length(data["I"]["Curr"],param["ny"]) # Store Length at time value
    
    # if we are at a frame we would like to plot, store the frame
    if telap%dtWindow == 0:
        print(f'Storing Frame at {telap}s')
        S[:,:,n] = data["S"]["Curr"]
        I[:,:,n] = data["I"]["Curr"]
        T[:,:,n] = data["T"]["Curr"]
        n+=1    # increment the frame counter
        
    f._global(data,param,dt) # Calculate the next time step
    
    t+=1            # Increment the length index
    telap=t*dt      # Increment the simulation time elapsed
```
The main function of this while loop is to calculate the values in the next time step given current time step. The if statement stores the images for plotting in the `S,I,T` variables of the `data` dictionary. The `f._global()` function call passes the `data, and param` dictionaries and the time step for the calculation (dt) into the `_global()` function of the `PDES.py` (defined as `f`) module where they are defined as **global variables**. This is where the calculation actually occur.
<hr>

**Note on Global Variables:** *I chose to define the variables globally (meaning they can be referenced and modified in each function) to eliminate the need to pass parameters into each function within the `PDES.py` file. This simulation can be done without doing so. For an example of this, please refer to the `functions.py` file within the `Summer 2020` folder.*
<hr>

## 3. Doing Some Math (PDES.py)
The main function of `PDE.py` file is to take the values  of our variables in our current time step and use them to calculate the values at the next time step.
```python
def pde():
    data["I"]["Curr"] = ice() # decide if something should be ice
    
    # Combine Diffusion and Advection on S
    [Fpxy,Fxx,Fyy,Fxy,Fyx] = diffusionbytemp("S")
    [Fxvel,Fyvel] = advectionbysalt("S",Fpxy)
    data["S"]["Curr"] = combine("S",Fxx,Fyy,Fxy,Fyx,Fxvel,Fyvel)
    
    # Combine Diffusion and Advection on T
    [Fpxy,Fxx,Fyy,Fxy,Fyx] = diffusion("T")
    [Fxvel,Fyvel] = advection("T",Fpxy)
    data["T"]["Curr"] = combine("T",Fxx,Fyy,Fxy,Fyx,Fxvel,Fyvel)
    
    # Set Constant Conditions
    spotwidth=int(param["nx"]*0.08) # This is half width in steps
    spotleft=int(np.round(param["nx"]/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(param["nx"]/2))+spotwidth  # Determine the right edge
    data["S"]["Curr"][spotleft:spotright,-1:] = param["brineSalinity"]
    data["T"]["Curr"][spotleft:spotright,-1:] = param["brineTemperature"]
```
This section of code can be broken down into three separate chunks:
<hr>

### 1. Decide if Something Should be Ice

The first step of our calculation is to determine if Ice is formed. This is done within the `ice()` function:

```python
def ice():
    # Load in Value of Arrays
    S = data["S"]["Curr"]
    I = data["I"]["Curr"]
    T = data["T"]["Curr"]
    
    I[:,:] = 0  # "Melt" all existing ice
 
    thresh = -1*10**(-10) # Temperature threshold
    result = np.where(T < thresh)
    x = result[0] # Gets a list of all the x indices
    y = result[1] # Gets a list of all the y indices
    ratio = param["ratio"] # Salinity(psu)/Temperature(◦C) ratio
    for i, n in enumerate(x):
        temp = T[x[i],y[i]]
        conc = S[x[i],y[i]]
        if (temp*ratio) < conc:
            I[x[i],y[i]] = 1
    
    return I
```
This function first loads in the values of `Salt, Ice, and Temperature` concentrations from the `data` dictionary defined in the `main.py` function. We then set all values of the `Ice` array to zero. This essentially "melts" all of the ice. We then find all values in the temperature array that are less than our `thresh` value. We then multiply the temperature at each of these points by the `ratio` value defined in the `param` dictionary within the `main.py` file.The current temperature multiplied by the `Salinity(psu)/Temperature(◦C)` ratio should give us the maximum concentration of salt we can have for freezing to still occur. If this ratio is less than the concentration value, we determine that ice has formed and change the index value within `I` to 1 to indicate Ice has formed. We leave index values as 0 if Ice is not formed.
<hr>

### 2. Combine Diffusion and Advection
<hr>
Directly after choosing if Ice should be formed, we begin performing Diffusion and Advection on both our temperature and salt concentrations. We do not do this for ice as ice does not diffuse since it is a solid and we define Ice as a boolean value meaning you are either ice or not ice and no value in between.

For the Temperature `diffusion` and `advection` we will use the following lines of function calls to calculate the change each variable has in both directions. 

```python
# Combine Diffusion and Advection on T
    [Fpxy,Fxx,Fyy,Fxy,Fyx] = diffusion("T")
    [Fxvel,Fyvel] = advection("T",Fpxy)
    data["T"]["Curr"] = combine("T",Fxx,Fyy,Fxy,Fyx,Fxvel,Fyvel)
```
**Diffusion:**

Diffusion is a movement of a concentration of something (in our use case: a fluid) from a region of higher concentration to a region of lower concentration:
<div align="center">
    <img src="../3D Models/Diffusion.gif" alt="drawing" width="300"/>
    <div> Figure 1: Python Simulation of Diffusion of a Square Fluid Concentration </div>
</div>

We implemented this in Python by first padding the array (generally referred to as `fp` ) with the second and second to last items such that there is no flux in the x direction (the last item has the same element on both sides).

```python
    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))
    fpxy = np.vstack((fpy[1,:], fpy, fpy[-2,:]))
```  
We then implemented a `finite differences method`:
<!-----------LATEX IN HTML----------->
<div align ="center"> 
    <img src="https://latex.codecogs.com/gif.latex?\dpi{125}&space;\bg_black&space;\fn_jvn&space;\boxed{\Delta[f](x)=f(x+1)-f(x).}"/>
</div>
<!-------------------------------------->

to define a partial differential equations for diffusion of a general variable `fp`:
<!-----------LATEX IN HTML----------->
<div align ="center"> 
    <img src="https://latex.codecogs.com/gif.latex?\dpi{125}&space;\bg_black&space;\fn_jvn&space;\boxed{fxx=dt*(df/res**2)*(fpxy[:,2:]+fpxy[:,0:-2]-2*fpxy[:,1:-1])}"/>
</div>
<!-------------------------------------->

```python
    # Perform finite differences. (Diffusion)
    fxx=dt*(df/res**2)*(fpxy[:,2:]+fpxy[:,0:-2]-2*fpxy[:,1:-1])
    fyy=dt*(df/res**2)*(fpxy[2:,:]+fpxy[0:-2,:]-2*fpxy[1:-1,:])
    fxx=fxx[1:-1,:] # Remove extra rows
    fyy=fyy[:,1:-1] # Remove extra rows
   
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15 # Set ff=1 to turn the fudge-factor off.
    
    # Diagonal terms
    fxy=dt*((df/(res**2+res**2))*(fpxy[2:,2:]+fpxy[0:-2,0:-2]-2*fpxy[1:-1,1:-1]))/ff
    fyx=dt*((df/(res**2+res**2))*(fpxy[0:-2,2:]+fpxy[2:,0:-2]-2*fpxy[1:-1,1:-1]))/ff
```

**Advection:**
<div align="center">
    <img src="../3D Models/AdvecDiff.gif" alt="drawing" width="300"/>
    <div> Figure 2: Python Simulation of Diffusion and Advection of a Square Fluid Concentration </div>
</div>


### 3. Set Constant Conditions
<hr>

We then set constant conditions for the simulation. We do this by defining a small square section at the top of our system where the brine is injected at a constant salt concentration and temperature defined within the  `brineSalinity` and `brineTemperature` values of the `param` dictionary of the `main.py` file:
```python
 # Set Constant Conditions
spotwidth=int(param["nx"]*0.08) # This is half width in steps
spotleft=int(np.round(param["nx"]/2))-spotwidth   # Determine the left edge
spotright=int(np.round(param["nx"]/2))+spotwidth  # Determine the right edge
data["S"]["Curr"][spotleft:spotright,-1:] = param["brineSalinity"]
data["T"]["Curr"][spotleft:spotright,-1:] = param["brineTemperature"]
```

## 4. Making Pretty Pictures (plot.py)
<hr>

By combining these methods, we can then create a "flip book" video of the different time steps of our calculation. These and videos are produced within the `Plots.py` file.
<div align="center">
    <img src="../3D Models/results.gif" alt="drawing" width="300"/>
    <div> Figure 3: Final output of the Python Simulation.</div>
</div>
As we can see, this results in a root t growth rate of the brinicle with a coefficient value of roughly 0.14.
<!-----------LATEX IN HTML----------->
<div align ="center"> 
    <img src="https://latex.codecogs.com/gif.latex?\dpi{125}&space;\bg_black&space;\fn_jvn&space;\boxed{L(t)\approx0.139\sqrt{t}}"/>
</div>
<!-------------------------------------->

## 5. Note on future work
- the data dictionary can be defined less sloppily.

- my implementation allows for more data about the variables to be stored but it's not needed currently.

- The Global variables could present some issues. Maybe remove this implementation.