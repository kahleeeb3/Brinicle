# -*- coding: utf-8 -*-
"""
Integrating ODEs using scipy

Created on Tue May 25 20:11:13 2021

@author: Caleb
"""
import numpy as np # Import numerical methods library
import matplotlib.pyplot as plt # Import plotting library
from scipy.integrate import odeint

#object falling
def f(Y,t,g):
    y,v = Y # Unpack vector Y = (y,v) of current state of the system
    dydt = v # Time rate of change of y-position
    dvdt = -g # Acceleration from Newton's second law
    dYdt = [dydt, dvdt] # Pack rates of change of dynamical variables into single vector
    return dYdt

y0 = 3 # Define initial position
v0 = 18 # Define initial velocity
g = 9.8 # Define gravitational acceleration
T = np.linspace(0,5) # Defines vector of times (1000 points between endpoints by default)

Y = odeint(f,[y0,v0],T,args=(g,)); # Integrate ODE numerically. If there were more parameters, they would go in the tuple with g.

# Plot results
plt.plot(T,Y[:,0],'.') # Plot position (stored in the first (0th) column of Y) vs. time
plt.xlabel('time [s]')
plt.ylabel('position [m]');