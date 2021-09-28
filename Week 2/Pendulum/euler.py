# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

%matplotlib inline
import numpy as np # Import numerical methods library
import matplotlib.pyplot as plt # Import plotting library


#  iterates one step forward in time using Euler's method.
def EulerStep(x0,deltat,f):
    return x0 + f(x0)*deltat

#  iterate EulerStep over and over until we hit a maximum time  𝑇
def Euler(x0,t0,T,deltat,f):
    # Set initial time and x value
    t = t0
    x = x0
    
    # Do Euler steps until t hits T
    while t < T:
        x = EulerStep(x,deltat,f) # Update x
        t = t + deltat # Update time
    
    return x

# defines the ODE used in our system
def Decay(x):
    y = -x
    return y

# Get the values
dt = [1,1e-1,1e-2,1e-3,1e-4,1e-5] # define the step size
calc = [0 for i in range(6)] # stores value returned by Euler function
for i, x in enumerate(dt):
    calc[i] = Euler(1,0,1,dt[i],Decay)

# Calculate the Error
exact = np.exp(-1) # exact solution to the ODE
error = [0 for i in range(6)]
for i, x in enumerate(calc):
    error[i] = abs(exact-calc[i])
    
# Plot
plt.scatter(np.log(dt),np.log(error))
plt.xlabel('Step Size')
plt.ylabel('Error')