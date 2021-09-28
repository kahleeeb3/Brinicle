#%% Imports
from PIL import Image, ImageOps
import numpy as np
import cv2
import os

# Imports for Output Function
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit

#%% PDE Functions
# pde function for no diffusion barrier
def diffusion(dt,res,xvel,yvel,df,fp):

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

# pde function for diffusion barrier
def diffusionb(dt,res,xvel,yvel,wp,df,fp):

    # Set up zero flux bc's in x.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.    
    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))
    fpxy = np.vstack((fpy[1,:], fpy, fpy[-2,:]))
    
    wpy = np.column_stack((wp[:,1],wp,wp[:,-2]))
    wpxy = np.vstack((wpy[1,:], wpy, wpy[-2,:]))
            
    # Perform finite differences. (Diffusion)
    fxx=dt*(1/res**2)*(df*(1-wpxy[:,1:-1]))*( (1-wpxy[:,2:])*(fpxy[:,2:]-fpxy[:,1:-1])+(1-wpxy[:,0:-2])*(fpxy[:,0:-2]-fpxy[:,1:-1]))
    fyy=dt*(1/res**2)*(df*(1-wpxy[1:-1,:]))*( (1-wpxy[2:,:])*(fpxy[2:,:]-fpxy[1:-1,:])+(1-wpxy[0:-2,:])*(fpxy[0:-2,:]-fpxy[1:-1,:]))
    fxx=fxx[1:-1,:] # Remove extra rows
    fyy=fyy[:,1:-1] # Remove extra rows
   
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15 # Set ff=1 to turn the fudge-factor off.
    
    # Diagonal terms
    fxy=dt*((1/(res**2+res**2))*(df*(1-wpxy[1:-1,1:-1]))*( (1-wpxy[2:,2:])*(fpxy[2:,2:]-fpxy[1:-1,1:-1])+(1-wpxy[0:-2,0:-2])*(fpxy[0:-2,0:-2]-fpxy[1:-1,1:-1])))/ff
    fyx=dt*((1/(res**2+res**2))*(df*(1-wpxy[1:-1,1:-1]))*( (1-wpxy[0:-2,2:])*(fpxy[0:-2,2:]-fpxy[1:-1,1:-1])+(1-wpxy[2:,0:-2])*(fpxy[2:,0:-2]-fpxy[1:-1,1:-1])))/ff
        
    return [fpxy,fxx,fyy,fxy,fyx]

# Diffusion Rate Depending on the Temperature
"""
This Doesn't Account for Ice Diffusion Barrier? It probably needs to.
"""
def diffusionbytemp(dt,res,xvel,yvel,df,fp,Tp):
    # https://dtrx.de/od/diff/
    # Diffusion rate of Water 1.05 to 0.187 μm^2/ms for water in temp of 0.35 °C to –30.65 °C
    # μm^2/ms = 10^(–3) mm^2/s
    # However, we assume constant rate of diffusion. Is this accpetable?
    # (D2-D1)/(T2-T1) = (0.187-1.05)/(-30.65-0.35) * 10**(-3)
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
    """
    Why are we doing 1 - diffRate?
    diffRate is around 0.000556774 - Very Small
    """
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
    
def advection(dt,res,xvel,fpxy,yvel,vel_f):
    # Caclulate the advection.
    #fyvel=dt*(vel_f*yvel/res)*(fpxy[:,0:-2]-fpxy[:,1:-1]) # Upward Advection
    fxvel=dt*(vel_f*xvel/res)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fyvel=dt*(vel_f*yvel/res)*(fpxy[:,2:]-fpxy[:,1:-1])
    fyvel=fyvel[1:-1,:] # Remove extra rows
    fxvel=fxvel[:,1:-1] # Remove extra rows
    
    return [fxvel,fyvel]

# pde function for diffusion barrier
"""
How Do I Do This?
"""
def advectionbysalt(dt,res,xvel,fpxy,yvel,vel_f,up):
    # Get the Advection Matrix from the concentration of salt at each point
    # With 32 PSU what should advection speed be?
    # With 60 PSU what should advection speed be?
    advRate = yvel * up
    zeroMatrix = 0*advRate
    advRate = np.column_stack((advRate[:,0],advRate,advRate[:,-1]))
    #advRate = np.vstack((advRate[1,:], advRate, advRate[-2,:]))
    # (above - at) - (at - below)
    # If less than 0, return 0 with np.maximum()
    advRatein = np.maximum(advRate[:,2:]-advRate[:,1:-1],zeroMatrix)
    advRateout = np.maximum(advRate[:,1:-1]-advRate[:,0:-2],zeroMatrix)
    #advRate = (advRatein-advRateout) / (60-32) # Unitless
    
    # fxvel=dt*(vel_f*xvel/res)*(fpxy[0:-2,:]-fpxy[1:-1,:])  
    # fyvel=dt*(vel_f*yvel/res)*(fpxy[:,2:]-fpxy[:,1:-1])
    fxvel= zeroMatrix #dt*(vel_f*xvel/res)*(fpxy[0:-2,:]-fpxy[1:-1,:])  
    fyvel=dt*(vel_f/res)*(advRatein*(fpxy[1:-1,2:])-advRateout*fpxy[1:-1,1:-1])
    
    
    #fyvel=fyvel[1:-1,:] # Remove extra rows
    #fxvel=fxvel[:,1:-1] # Remove extra rows
    
    return [fxvel,fyvel/(60-32)]

def combine(fp,fxx,fyy,fxy,fyx,fxvel,fyvel):      
    # Combine diffusion with advection and kinetics using the forward Euler algorithm.
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
        
    return fp

# Calls pdef1() and pdef2() above
def pde(data,diff,vel,xvel,yvel,nx,ny,dt,res,telap,brineSalinity,brineTemperature,ratio):    
    
    #load the data
    u = data[0][:,:]
    T = data[1][:,:]
    w = data[2][:,:]
    
    # decide if something should be ice
    w[:,:] = ice(u,T,w,ratio)
    
    # Pass each array into their respective function
    
    # Combine Diffusion and Advection on u
    [upxy,uxx,uyy,uxy,uyx] = diffusionbytemp(dt,res,xvel,yvel,diff[0],u,T)
    #[uxvel,uyvel] = advection(dt,res,xvel,upxy,yvel,vel[0])
    [uxvel,uyvel] = advectionbysalt(dt,res,xvel,upxy,4.655,vel[0],u)
    u = combine(u,uxx,uyy,uxy,uyx,uxvel,uyvel)
    
    # Combine Diffusion and Advection on T
    [Tpxy,Txx,Tyy,Txy,Tyx] = diffusion(dt,res,xvel,yvel,diff[1],T)
    [Txvel,Tyvel] = advection(dt,res,xvel,Tpxy,yvel,vel[1])
    T = combine(T,Txx,Tyy,Txy,Tyx,Txvel,Tyvel)
    
    # Set Constant Conditions
    spotwidth=int(nx*0.08) # This is half width in steps
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    u[spotleft:spotright,-1:] = brineSalinity
    T[spotleft:spotright,-1:] = brineTemperature
        
    # Return the resulting array
    return [u, T, w]

#%% Ice Calculations

# Decides if there is ice or not
def ice(u,T,w,ratio):
    
    """
    The current temperature multiplied by the 
    Salinity(psu)/Temperature(◦C) ratio should
    give us the maximum conentration of salt 
    we can have for freezing to still occur
    """
    
    #w[:,:] = 0  # "Melt" all existing ice
        
    # Didnt set value to < 0 b/c there is a 
    thresh = -1*10**(-10) #accounts for values that are basically zero but technically are less than
    result = np.where(T < thresh) # Finds Where the Temperature is Less than 0
    x = result[0] # Gets a list of all the x indices
    y = result[1] # Gets a list of all the y indices
    for i, n in enumerate(x):
        temp = T[x[i],y[i]]
        conc = u[x[i],y[i]]
        if (temp*ratio) < conc:
            w[x[i],y[i]] = 1
    
    return w

def length(lp,wp,ny):
    w_length= np.sum(wp[:,:],axis=0)
    w_nonzeros=np.flatnonzero(w_length)
    try:
        w_end=ny-w_nonzeros[0]
    except:
        w_end=0
    return w_end/100
    
#%% Plots

def func(x, a):
    return a * np.sqrt(x)
        
# The new plot function (unstable)
def plot(folder,tEnd,Lx, Ly, res, data, l, dt, dtWindow, time):
    
    # get data from matrix
    x = np.arange(0,Lx,res) # Create the x data for plotting
    y = np.arange(0,Ly,res) # Create the y data for plotting
    t = math.floor(time/dt)
    Y, X = np.meshgrid(y,x) # Create the X,Y matrices for plotting
    
    # define plots
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 3), (0, 0))
    ax2 = plt.subplot2grid((2, 3), (0, 1))
    ax3 = plt.subplot2grid((2, 3), (0, 2))
    ax4 = plt.subplot2grid((2, 1), (1, 0))#, colspan=2)
    
    # aspect ratio
    ax1.set_aspect('equal')
    ax2.set_aspect('equal')
    ax3.set_aspect('equal')
    #ax4.set_aspect(1)
    
    # plot functions
    cf1 = ax1.contourf(X,Y,data[0],cmap=cm.coolwarm)
    cf2 = ax2.contourf(X,Y,data[1],cmap=cm.coolwarm)
    cf3 = ax3.contourf(X,Y,data[2],cmap=cm.coolwarm)
    xdata = l[0,:t]
    ydata = l[1,:t]
    ax4.plot(xdata,ydata, label='Brinicle Length') # Length Plot
    popt, pcov = curve_fit(func, xdata, ydata) # Perform Curve Fit
    ax4.plot(xdata, func(xdata, *popt), '--',label='f(x)=a√(x), a=%5.3f' % tuple(popt))
    ax4.legend()
    #ax4.plot(l[0,0:t],np.sqrt(l[0,0:t])) # root(t) plot
    
    fig.colorbar(cf1, ax=ax1) # Add the colorbar
    fig.colorbar(cf2, ax=ax2) # Add the colorbar
    fig.colorbar(cf3, ax=ax3) # Add the colorbar
    
    
    size = 9
    
    ax1.set_xlabel('x (mm)') # Label the x axis
    ax1.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration\n of Salt at {0:0.1f}s'.format(time)
    ax1.set_title(utitle,fontsize=size)
    
    ax2.set_xlabel('x (mm)') # Label the x axis
    ax2.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration\n of Temp at {0:0.1f}s'.format(time)
    ax2.set_title(utitle,fontsize=size)
    
    ax3.set_xlabel('x (mm)') # Label the x axis
    ax3.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration\n of Ice at {0:0.1f}s'.format(time)
    ax3.set_title(utitle,fontsize=size)
    
    ax4.set_xlabel('x (s)') # Label the x axis
    ax4.set_ylabel('y (cm)') # Label the y axis
    utitle = 'Length of Ice at {0:0.1f}s'.format(time)
    ax4.set_title(utitle,fontsize=size)
    
    # text display
    """
    eqtext = "du/dt = auv\ndv/dt = -buv+cuv\ndw/dt = duv"
    vartext = "u - salt water \nv - brine \nw - ice \n"
    text = vartext + eqtext
    textx, texty, size = 1.05, 0, 8
    plt.text(textx, texty, text, fontsize=size)
    """
    
    #plt.show()
    plt.tight_layout()
    fname = 'images/Brinicle_%06d.png' % t # Create the file name for each plot
    print('Saving frame', fname) # Print the status update
    fig.savefig(fname, dpi=300) # Save the image
    plt.close(fig) # Close the image so it doesn't show while the code is running
    plt.clf() # Clear the figrue to save memory    

#%% Image Handling

# Returns an array from an image file named image.png
def image(name,x, y, angle):
    
    im = Image.open(name)    #Open the image
    
    im = im.resize((y, x))  #Change Image Size
    im = im.rotate(angle)   #Rotate the Image
    im = im.convert('L')    #Convert to Black & White
    im = ImageOps.invert(im)   #Invert the image
    
    
    #im.save('image2.png') #Save the new image
    
    #Return the image as an array
    data = np.asarray(im)
    return data

#Converts image folder to a video
def export(image_folder):
    
    fps = 10
    print(f'Exporting: {image_folder} to {fps}fps Video')
    video_name = f'{image_folder}.avi'
    

    
    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape
    
    video = cv2.VideoWriter(video_name, 0, fps, (width,height))
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))
    
    cv2.destroyAllWindows()
    video.release()
    print(f'Completed: Saved to {image_folder}.avi')
    