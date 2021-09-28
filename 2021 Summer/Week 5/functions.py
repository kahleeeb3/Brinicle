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
    
#%% Length Calculations

def length(lp,wp,res):
    #loop through the indecies of the middle
    # w is size (60,170,t)
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
 
#%% Plots
        
# The new plot function (unstable)
def plot(folder,tEnd,Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, time):
    
    # get data from matrix
    x = np.arange(0,Lx,dx) # Create the x data for plotting
    y = np.arange(0,Ly,dy) # Create the y data for plotting
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
    utitle = 'Concentration\n of u at {0:0.1f}s'.format(time)
    ax1.set_title(utitle,fontsize=size)
    
    ax2.set_xlabel('x (mm)') # Label the x axis
    ax2.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration\n of v at {0:0.1f}s'.format(time)
    ax2.set_title(utitle,fontsize=size)
    
    ax3.set_xlabel('x (mm)') # Label the x axis
    ax3.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration\n of w at {0:0.1f}s'.format(time)
    ax3.set_title(utitle,fontsize=size)
    
    ax4.set_xlabel('x (s)') # Label the x axis
    ax4.set_ylabel('y (cm)') # Label the y axis
    utitle = 'Length of w at {0:0.1f}s'.format(time)
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

# pde function for no diffusion barrier
def pdef(dt,dx,dy,xvel,yvel,vel_f,df,fp):

    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))
    fpxy = np.vstack((fpy[1,:], fpy, fpy[-2,:]))
            
    fxx=dt*(df/dx**2)*(fpxy[:,2:]+fpxy[:,0:-2]-2*fpxy[:,1:-1])
    fxx=fxx[1:-1,:]
    fyy=dt*(df/dy**2)*(fpxy[2:,:]+fpxy[0:-2,:]-2*fpxy[1:-1,:])
    fyy=fyy[:,1:-1]
   
    ff=15
    
    fxy=dt*((df/(dx**2+dy**2))*(fpxy[2:,2:]+fpxy[0:-2,0:-2]-2*fpxy[1:-1,1:-1]))/ff
    fyx=dt*((df/(dx**2+dy**2))*(fpxy[0:-2,2:]+fpxy[2:,0:-2]-2*fpxy[1:-1,1:-1]))/ff
    
    fxvel=dt*(vel_f*xvel/dx)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fxvel=fxvel[:,1:-1]
    
    #fyvel=dt*(vel_f*yvel/dy)*(fpxy[:,0:-2]-fpxy[:,1:-1]) # Upward Advection
    fyvel=dt*(vel_f*yvel/dy)*(fpxy[:,2:]-fpxy[:,1:-1])
    fyvel=fyvel[1:-1,:]
            
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
        
    return fp

# Calls pdef1() and pdef2() above
def pde(data,diff,vel,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d):
    
    # Pass each array into their respective function
    u = pdef(dt,dx,dy,xvel,yvel,vel[0],diff[0],data[0])  
    v = pdef(dt,dx,dy,xvel,yvel,vel[1],diff[1],data[1])
    w = pdef(dt,dx,dy,xvel,yvel,vel[2],diff[2],data[2])
    
    # Apply the PDE to the new Arrays
    u = u + dt*(-a*u*v)
    v = v + dt*(-b*u*v+c*u*v)
    w = w + dt*(d*u*v)
    
    # Constant Conditions
    spotwidth=8 # This is half width in steps
    spotheight =1
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    v[spotleft:spotright,ny-spotheight:ny]= 1      # Constant Brine Value
        
    # Return the resulting array
    return [u, v, w]