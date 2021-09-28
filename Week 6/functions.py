#%% Imports
from PIL import Image, ImageOps
import numpy as np
import cv2
import os

# Imports for Output Function
import math
import matplotlib.pyplot as plt
from matplotlib import cm

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
    
#%% Ice Calculations

# Decides if there is ice or not
def ice(u,T,w,ratio):
    
    """
    The current temperature multiplied by the 
    Salinity(psu)/Temperature(◦C) ratio should
    give us the maximum conentration of salt 
    we can have for freezing to still occur
    """
    
    w[:,:] = 0  # "Melt" all existing ice
        
    # Didnt set value to < 0 b/c there is a lot of values that are basically zero but technically are less than
    result = np.where(T < -0.00001) # Finds Where the Temperature is Less than 0
    x = result[0] # Gets a list of all the x indices
    y = result[1] # Gets a list of all the y indices
    for i, n in enumerate(x):
        temp = T[x[i],y[i]]
        conc = u[x[i],y[i]]
        if (temp*ratio) < conc:
            w[x[i],y[i]] = 1
    
    """
    w[:,:] = 0  #Remove all existing Ice
    temp = -2              # defines the temperature for ice formation
    salt_conc = 5          # Freezing temp of water as a function of salinity
    
    result = np.where(u > salt_conc) # Finds places in the array where salt concentration is a certain value
    x = result[0] # Gets a list of all the x indices
    y = result[1] # Gets a list of all the y indices
    for i, n in enumerate(x):
        if T[x[i],y[i]] < temp:
            w[x[i],y[i]] = 1
    """
    
    return w

def length(lp,wp,ny):
    w_length= np.sum(wp[:,:],axis=0)
    w_nonzeros=np.flatnonzero(w_length)
    try:
        w_end=ny-w_nonzeros[0]
    except:
        w_end=0
    return w_end/10
    """
    l = 0
    w = len(wp[:,0]) # find how wide the array is
    h = len(wp[0,:]) # find how tall the array is
    c = math.floor(w/2) # finds the center index
    
    #loop through every row in the middle column
    for x in range(h):
        if wp[c,x] > 0.01:
            l = x
            break
        if x == h-1:
            l = h # no ice was found
            
    l = (h - l) #*res # convert index into measurement
    return l/10
    """
 
#%% Plots
        
# The new plot function (unstable)
def plot(folder,tEnd,Lx, Ly, res, u, v, w, l, dt, dtWindow, time):
    
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
    cf1 = ax1.contourf(X,Y,u,cmap=cm.coolwarm)
    cf2 = ax2.contourf(X,Y,v,cmap=cm.coolwarm)
    cf3 = ax3.contourf(X,Y,w,cmap=cm.coolwarm)
    ax4.plot(l[0,:t],l[1,:t]) # Length Plot
    ax4.plot(l[0,:t],np.sqrt(l[0,:t])) # root(t) plot
    
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

#%% PDE Functions

# pde function for diffusion barrier
def pdedb(dt,res,xvel,yvel,wp,vel_f,df,fp):

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
    
    # Caclulate the advection.
    #fyvel=dt*(vel_f*yvel/res)*(fpxy[:,0:-2]-fpxy[:,1:-1]) # Upward Advection
    fxvel=dt*(vel_f*xvel/res)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fyvel=dt*(vel_f*yvel/res)*(fpxy[:,2:]-fpxy[:,1:-1])
    fyvel=fyvel[1:-1,:] # Remove extra rows
    fxvel=fxvel[:,1:-1] # Remove extra rows
            
    # Combine diffusion with advection and kinetics using the forward Euler algorithm.
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
        
    return fp

# pde function for no diffusion barrier
def pded(dt,res,xvel,yvel,vel_f,df,fp):

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
    
    # Caclulate the advection.
    #fyvel=dt*(vel_f*yvel/res)*(fpxy[:,0:-2]-fpxy[:,1:-1]) # Upward Advection
    fxvel=dt*(vel_f*xvel/res)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fyvel=dt*(vel_f*yvel/res)*(fpxy[:,2:]-fpxy[:,1:-1])
    fyvel=fyvel[1:-1,:] # Remove extra rows
    fxvel=fxvel[:,1:-1] # Remove extra rows
            
    # Combine diffusion with advection and kinetics using the forward Euler algorithm.
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
        
    return fp

# Calls pdef1() and pdef2() above
def pde(data,diff,vel,xvel,yvel,nx,ny,dt,res,telap,brineSalinity,brineTemperature,ratio):    
    
    #load the data
    u = data[0]
    T = data[1]
    w = data[2]
    
    # decide if something should be ice
    w[:,:] = ice(u,T,w,ratio)
    
    # Pass each array into their respective function
    u = pdedb(dt,res,xvel,yvel,w,vel[0],diff[0],u)  
    T = pded(dt,res,xvel,yvel,vel[1],diff[1],T)
    
    # Shift Rows Down
    spotwidth=8 # This is half width in steps
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    
    t = math.floor(telap/dt)
    advel = math.floor((res/dt)/yvel)
    if t%advel == 0:
        
        u[spotleft:spotright,0:-1] = u[spotleft:spotright,1:]
        T[spotleft:spotright,0:-1] = T[spotleft:spotright,1:]
    
        # Set Constant Conditions
        u[spotleft:spotright,-1:] = brineSalinity
        T[spotleft:spotright,-1:] = brineTemperature
        
    # Return the resulting array
    return [u, T, w]