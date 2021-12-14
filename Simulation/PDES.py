#%% Imports
#from PIL import Image, ImageOps
import numpy as np
#import cv2
#import os

# Imports for Output Function
#import math
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from scipy.optimize import curve_fit

#%% Define Global Variables
def _global(data_p,param_p,dt_p):
    # create global variables
    global data, param
    data = data_p
    param = param_p
    
    # CONSTANTS
    global dt,xvel,yvel,res
    res = param["res"]
    xvel = param["xvel"]
    yvel = param["yvel"]
    dt = dt_p
    
    pde() # Run the PDE Function
    
    return data

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
     
def ice():
    """
    The current temperature multiplied by the Salinity(psu)/Temperature(◦C) ratio should
    give us the maximum conentration of salt we can have for freezing to still occur
    """
    # Load in Vlaue of Arrays
    S = data["S"]["Curr"]
    I = data["I"]["Curr"]
    T = data["T"]["Curr"]
    
    I[:,:] = 0  # "Melt" all existing ice
        
    # Didnt set value to < 0 b/c there is a 
    thresh = -1*10**(-10) #accounts for values that are basically zero but technically are less than
    result = np.where(T < thresh) # Finds Where the Temperature is Less than 0
    x = result[0] # Gets a list of all the x indices
    y = result[1] # Gets a list of all the y indices
    ratio = param["ratio"] # Salinity(psu)/Temperature(◦C) ratio
    for i, n in enumerate(x):
        temp = T[x[i],y[i]]
        conc = S[x[i],y[i]]
        if (temp*ratio) < conc:
            I[x[i],y[i]] = 1
    
    return I

def length(wp,ny):
    w_length= np.sum(wp[:,:],axis=0)
    w_nonzeros=np.flatnonzero(w_length)

    try:
        w_end=ny-w_nonzeros[0]
    except:
        w_end=0
    
    return w_end/100

def combine(var_name,fxx,fyy,fxy,fyx,fxvel,fyvel):      
    # Combine diffusion with advection and kinetics using the forward Euler algorithm.
    fp = data[var_name]["Curr"]
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
        
    return fp

#%% Diffusion
def diffusion(var_name):
    fp = data[var_name]["Curr"]
    df = data[var_name]["diff"]
    
    # Set up zero flux bc's in x.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.    
    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))
    fpxy = np.vstack((fpy[1,:], fpy, fpy[-2,:]))
            
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
    
    return [fpxy,fxx,fyy,fxy,fyx]

     #diffusionbytemp(dt,res,xvel,yvel,diff[0],u,T)
# def diffusionbytemp(dt,res,xvel,yvel,df,fp,Tp):
def diffusionbytemp(var_name):
    """
    https://dtrx.de/od/diff/
    - Diffusion rate of Water 1.05 to 0.187 μm^2/ms for water in temp of 0.35 °C to –30.65 °C
    - μm^2/ms = 10^(–3) mm^2/s
    - However, we assume constant rate of diffusion. Is this accpetable?
    - (D2-D1)/(T2-T1) = (0.187-1.05)/(-30.65-0.35) * 10**(-3)
    """
    Tp = data["T"]["Curr"] # Temperature Value
    fp = data[var_name]["Curr"] # Variable to Diffuse
    df = data[var_name]["diff"] # Diffusion Rate of Variable
    
    
    rate = (0.187-1.05)/(30.65+0.35) * 10**(-3)
    diffRate = rate * Tp # This is the diffusion Array
    
    # Set up zero flux bc's in x.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.    
    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))
    fpxy = np.vstack((fpy[1,:], fpy, fpy[-2,:]))
    
    diffRatey = np.column_stack((diffRate[:,1],diffRate,diffRate[:,-2]))
    diffRatexy = np.vstack((diffRatey[1,:], diffRatey, diffRatey[-2,:]))
    
    #wpy = np.column_stack((wp[:,1],wp,wp[:,-2]))
    #wpxy = np.vstack((wpy[1,:], wpy, wpy[-2,:]))
            
    # Perform finite differences. (Diffusion)
    fxx=dt*(1/res**2)*(df*(1-diffRatexy[:,1:-1]))*( (1-diffRatexy[:,2:])*(fpxy[:,2:]-fpxy[:,1:-1])+(1-diffRatexy[:,0:-2])*(fpxy[:,0:-2]-fpxy[:,1:-1]))
    fyy=dt*(1/res**2)*(df*(1-diffRatexy[1:-1,:]))*( (1-diffRatexy[2:,:])*(fpxy[2:,:]-fpxy[1:-1,:])+(1-diffRatexy[0:-2,:])*(fpxy[0:-2,:]-fpxy[1:-1,:]))
    fxx=fxx[1:-1,:] # Remove extra rows
    fyy=fyy[:,1:-1] # Remove extra rows
   
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15 # Set ff=1 to turn the fudge-factor off.
    
    # Diagonal terms
    fxy=dt*((1/(res**2+res**2))*(df*(1-diffRatexy[1:-1,1:-1]))*( (1-diffRatexy[2:,2:])*(fpxy[2:,2:]-fpxy[1:-1,1:-1])+(1-diffRatexy[0:-2,0:-2])*(fpxy[0:-2,0:-2]-fpxy[1:-1,1:-1])))/ff
    fyx=dt*((1/(res**2+res**2))*(df*(1-diffRatexy[1:-1,1:-1]))*( (1-diffRatexy[0:-2,2:])*(fpxy[0:-2,2:]-fpxy[1:-1,1:-1])+(1-diffRatexy[2:,0:-2])*(fpxy[2:,0:-2]-fpxy[1:-1,1:-1])))/ff
        
    return [fpxy,fxx,fyy,fxy,fyx]

#%% advection
def advection(var_name,fpxy):
    vel_f = data[var_name]["vel"]
    
    # Caclulate the advection.
    #fyvel=dt*(vel_f*yvel/res)*(fpxy[:,0:-2]-fpxy[:,1:-1]) # Upward Advection
    fxvel=dt*(vel_f*xvel/res)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fyvel=dt*(vel_f*yvel/res)*(fpxy[:,2:]-fpxy[:,1:-1])
    fyvel=fyvel[1:-1,:] # Remove extra rows
    fxvel=fxvel[:,1:-1] # Remove extra rows
    
    return [fxvel,fyvel]

def advectionbysalt(var_name,fpxy):
    vel_f = data[var_name]["vel"]
    #advRate = yvel * data[var_name]["Curr"]
    advRate = 4.655 * data[var_name]["Curr"] # I dont remember why we did this?
    zeroMatrix = 0*advRate
    advRate = np.column_stack((advRate[:,0],advRate,advRate[:,-1]))
    
    advRatein = np.maximum(advRate[:,2:]-advRate[:,1:-1],zeroMatrix)
    advRateout = np.maximum(advRate[:,1:-1]-advRate[:,0:-2],zeroMatrix)
    
    fxvel= zeroMatrix 
    fyvel=dt*(vel_f/res)*(advRatein*(fpxy[1:-1,2:])-advRateout*fpxy[1:-1,1:-1])
    
    return [fxvel,fyvel/(60-32)] # Why did we do 60-32?

