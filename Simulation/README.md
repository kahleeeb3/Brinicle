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
    <img src="https://latex.codecogs.com/gif.latex?\dpi{75}&space;\bg_black&space;\fn_jvn&space;\boxed{x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}.}"/>
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
    <img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;\bg_black&space;\fn_jvn&space;\boxed{x_n = x_{n+1}.}"/>
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

```python
def pde():
```

## 4. Making Pretty Pictures (plot.py)