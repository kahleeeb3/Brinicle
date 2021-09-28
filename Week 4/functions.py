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
        
# The original plot function (stable)
def plotOG(folder,tEnd,Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, time):
    t = math.floor(time/dt)
    telap = time
    
    #Make the plots
    x = np.arange(0,Lx,dx) # Create the x data for plotting
    y = np.arange(0,Ly,dy) # Create the y data for plotting
    Y, X = np.meshgrid(y,x) # Create the X,Y matrices for plotting
    Z1 = u[:,:,t] # Pull the Z data matrix for plotting
    Z2 = v[:,:,t] # Pull the Z data matrix for plotting
    Z3 = w[:,:,t] # Pull the Z data matrix for plotting
    fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2,2) # Create the figure with subplots
    
    #cf1
    # Create the filled countour plot with colormap and manual levels
    #cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm,levels=np.arange(0,windowmax,res))
    cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm)
    fig.colorbar(cf1, ax=ax1) # Add the colorbar
    ax1.set_xlabel('x (mm)') # Label the x axis
    ax1.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Concentration of u \n at {0:0.1f}s'.format(telap)
    ax1.set_title(utitle) # Title the plot
    ax1.set_aspect('equal') # Make the aspect ratio equal
    
    #cf2
    # Create the filled countour plot with colormap and manual levels
    #cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm,levels=np.arange(0,windowmax,res))
    cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm)
    fig.colorbar(cf2, ax=ax2) # Add the colorbar
    ax2.set_xlabel('x (mm)') # Label the x axis
    ax2.set_ylabel('y (mm)') # Label the y axis
    vtitle = 'Concentration of v \n at {0:0.1f}s'.format(telap)
    ax2.set_title(vtitle) # Title the plot
    ax2.set_aspect('equal') # Make the aspect ratio equal
    
    #cf3
    # Create the filled countour plot with colormap and manual levels
    #cf3 = ax3.contourf(X,Y,Z3,cmap=cm.coolwarm,levels=np.arange(0,windowmax,res))
    cf3 = ax3.contourf(X,Y,Z3,cmap=cm.coolwarm)
    fig.colorbar(cf3, ax=ax3) # Add the colorbar
    ax3.set_xlabel('x (mm)') # Label the x axis
    ax3.set_ylabel('y (mm)') # Label the y axis
    wtitle = 'Concentration of w \n at {0:0.1f}s'.format(telap)
    ax3.set_title(wtitle) # Title the plot
    ax3.set_aspect('equal') # Make the aspect ratio equal
   
    #cf4
    # Create the filled countour plot with colormap and manual levels
    #cf4 = ax4.contourf(X,Y,Z4,cmap=cm.coolwarm,levels=np.arange(0,windowmax,res))
    #cf4 = ax4.contourf(X,Y,Z4,cmap=cm.coolwarm)
    ax4.plot(l[0,:t],l[1,:t])
    #ax4.plot(l[0,:t],np.sqrt(l[0,:t])) # root(t) plot
    ax4.set_xlabel('x (s)') # Label the x axis
    ax4.set_ylabel('y (cm)') # Label the y axis
    utitle = 'Length of w \n at {0:0.1f}s'.format(telap)
    ax4.set_title(utitle) # Title the plot
    
    #plt.text(textx, texty, text, fontsize=size)
    plt.tight_layout()

    #plt.show() # This shows the plots as the code is running
    fname = 'images/Brinicle_%06d.png' % t # Create the file name for each plot
    print('Saving frame', fname) # Print the status update
    fig.savefig(fname, dpi=300) # Save the image
    plt.close(fig) # Close the image so it doesn't show while the code is running
    plt.clf() # Clear the figrue to save memory
    
# The new plot function (unstable)
def plotnew(folder,tEnd,Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, time):
    
    # get data from matrix
    x = np.arange(0,Lx,dx) # Create the x data for plotting
    y = np.arange(0,Ly,dy) # Create the y data for plotting
    t = math.floor(time/dt)
    Y, X = np.meshgrid(y,x) # Create the X,Y matrices for plotting
    Z1 = u[:,:,t] # Pull the Z data matrix for plotting
    Z2 = v[:,:,t] # Pull the Z data matrix for plotting
    Z3 = w[:,:,t] # Pull the Z data matrix for plotting
    
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
    cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm)
    cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm)
    cf3 = ax3.contourf(X,Y,Z3,cmap=cm.coolwarm)
    ax4.plot(l[0,:t],l[1,:t])
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

# runs the functions above
def plot(folder,tEnd,Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, time):
    if time == 0:
         t = 0        # This sets the first time point for the calculation.
         telap = 0    # This sets the time elapsed for the simulation.
         while telap <= tEnd:  
             #plotOG(folder,tEnd, Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, telap)
             plotnew(folder,tEnd, Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, telap)
             t=t+math.floor(dtWindow/dt)            # Increment the storage counter
             telap=t*dt         # Increment the simulation time elapsed
    else:
        #plotOG(folder,tEnd, Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, time)
        plotnew(folder,tEnd, Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow, time)

#%% Old PDE Function

def ard_pde(up,vp,wp,du,dv,dw,vel_u,vel_v,vel_w,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d):
# This numerically integrates using finite differences and the forward
# Euler method.
        
    # Set up zero flux bc's in x.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
    upy = np.column_stack((up[:,1],up,up[:,-2]))
    vpy = np.column_stack((vp[:,1],vp,vp[:,-2]))
    wpy = np.column_stack((wp[:,1],wp,wp[:,-2]))
    # Make boundaries zero flux in y.
    upxy = np.vstack((upy[1,:], upy, upy[-2,:]))
    vpxy = np.vstack((vpy[1,:], vpy, vpy[-2,:]))
    wpxy = np.vstack((wpy[1,:], wpy, wpy[-2,:]))
            
    # Perform finite differences. (Diffusion)
    # Calculate differences
    uxx=dt*(du/dx**2)*(upxy[:,2:]+upxy[:,0:-2]-2*upxy[:,1:-1])
    vxx=dt*(dv/dx**2)*(vpxy[:,2:]+vpxy[:,0:-2]-2*vpxy[:,1:-1])
    wxx=dt*(dw/dx**2)*(wpxy[:,2:]+wpxy[:,0:-2]-2*wpxy[:,1:-1])
    
    uyy=dt*(du/dy**2)*(upxy[2:,:]+upxy[0:-2,:]-2*upxy[1:-1,:])
    vyy=dt*(dv/dy**2)*(vpxy[2:,:]+vpxy[0:-2,:]-2*vpxy[1:-1,:])
    wyy=dt*(dw/dy**2)*(wpxy[2:,:]+wpxy[0:-2,:]-2*wpxy[1:-1,:])
    
    uxx=uxx[1:-1,:] # Remove extra rows
    vxx=vxx[1:-1,:] # Remove extra rows
    wxx=wxx[1:-1,:] # Remove extra rows
    
    uyy=uyy[:,1:-1] # Remove extra columns    
    vyy=vyy[:,1:-1] # Remove extra columns  
    wyy=wyy[:,1:-1] # Remove extra columns
   
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15;  # Set ff=1 to turn the fudge-factor off.
    
    # Diagonal terms for u
    # Calculate differences
    uxy=dt*((du/(dx**2+dy**2))*(upxy[2:,2:]+upxy[0:-2,0:-2]-2*upxy[1:-1,1:-1]))/ff
    vxy=dt*((dv/(dx**2+dy**2))*(vpxy[2:,2:]+vpxy[0:-2,0:-2]-2*vpxy[1:-1,1:-1]))/ff
    wxy=dt*((dw/(dx**2+dy**2))*(wpxy[2:,2:]+wpxy[0:-2,0:-2]-2*wpxy[1:-1,1:-1]))/ff
    
    uyx=dt*((du/(dx**2+dy**2))*(upxy[0:-2,2:]+upxy[2:,0:-2]-2*upxy[1:-1,1:-1]))/ff
    vyx=dt*((dv/(dx**2+dy**2))*(vpxy[0:-2,2:]+vpxy[2:,0:-2]-2*vpxy[1:-1,1:-1]))/ff
    wyx=dt*((dw/(dx**2+dy**2))*(wpxy[0:-2,2:]+wpxy[2:,0:-2]-2*wpxy[1:-1,1:-1]))/ff
    
    # Caclulate the advection along the x axis.
    uxvel=dt*(vel_u*xvel/dx)*(upxy[0:-2,:]-upxy[1:-1,:])
    vxvel=dt*(vel_v*xvel/dx)*(vpxy[0:-2,:]-vpxy[1:-1,:])
    wxvel=dt*(vel_w*xvel/dx)*(wpxy[0:-2,:]-wpxy[1:-1,:])
    
    uxvel=uxvel[:,1:-1] # Remove extra rows
    vxvel=vxvel[:,1:-1] # Remove extra rows
    wxvel=wxvel[:,1:-1] # Remove extra rows
    
    # Caclulate the advection along the y axis.
    #uyvel=dt*(vel_u*yvel/dy)*(upxy[:,0:-2]-upxy[:,1:-1])
    #wyvel=dt*(vel_w*yvel/dy)*(wpxy[:,0:-2]-wpxy[:,1:-1])
    uyvel=dt*(vel_u*yvel/dy)*(upxy[:,2:]-upxy[:,1:-1])
    vyvel=dt*(vel_v*yvel/dy)*(vpxy[:,2:]-vpxy[:,1:-1])
    wyvel=dt*(vel_w*yvel/dy)*(wpxy[:,2:]-wpxy[:,1:-1])
    
    
    uyvel=uyvel[1:-1,:] # Remove extra rows
    vyvel=vyvel[1:-1,:] # Remove extra rows
    wyvel=wyvel[1:-1,:] # Remove extra rows
            
    # Combine diffusion with advection and kinetics using
    # the forward Euler algorithm.
    up=up+(uxx+uyy+uxy+uyx)+(uxvel+uyvel) + dt*(-a*up*vp)
    vp=vp+(vxx+vyy+vxy+vyx)+(vxvel+vyvel) + dt*(-b*up*vp+c*up*vp)
    wp=wp+(wxx+wyy+wxy+wyx)+(wxvel+wyvel) + dt*(d*up*vp)
    
    # Constant Conditions
    # Constant Brine Source in Center
    spotwidth=8 # This is half width in steps
    spotheight =1
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    vp[spotleft:spotright,ny-spotheight:ny]= 1      # Constant Brine Value
        
    # Return the resulting array
    return [up, vp, wp]

#%% New PDE Functions

# pde function that alters code for a diffusion barrier
def pdef_db(wp,dt,dx,dy,xvel,yvel,vel_f,df,fp):
    # This numerically integrates using finite differences and the forward Euler method.
    # Rates of diffusion matrix
    dfm = df/(1+wp)
    dfmy = np.column_stack((dfm[:,1],dfm,dfm[:,-2]))
    dfmxy = np.vstack((dfmy[1,:],dfmy, dfmy[-2,:]))
    
    # zero flux boundary conditions
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))
    fpxy = np.vstack((fpy[1,:],fpy,fpy[-2,:]))
    
    # Diffusion
    fxx=(dt/dx**2)*(dfmxy[:,2:]*fpxy[:,2:]+dfmxy[:,0:-2]*fpxy[:,0:-2]-2*dfmxy[:,1:-1]*fpxy[:,1:-1])  
    fxx=fxx[1:-1,:] # Remove extra rows
    fyy=(dt/dy**2)*(dfmxy[2:,:]*fpxy[2:,:]+dfmxy[0:-2,:]*fpxy[0:-2,:]-2*dfmxy[1:-1,:]*fpxy[1:-1,:])
    fyy=fyy[:,1:-1] # Remove extra columns
    
    # Fudge Factor
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15;  # Set ff=1 to turn the fudge-factor off.
    
    # Diffusion (diagonal terms)
    fxy=((dt/(dx**2+dy**2))*(dfmxy[2:,2:]*fpxy[2:,2:]+dfmxy[0:-2,0:-2]*fpxy[0:-2,0:-2]-2*dfmxy[1:-1,1:-1]*fpxy[1:-1,1:-1]))/ff
    fyx=((dt/(dx**2+dy**2))*(dfmxy[0:-2,2:]*fpxy[0:-2,2:]+dfmxy[2:,0:-2]*fpxy[2:,0:-2]-2*dfmxy[1:-1,1:-1]*fpxy[1:-1,1:-1]))/ff
    
    # Advection
    # x-axis
    fxvel=dt*(vel_f*xvel/dx)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fxvel=fxvel[:,1:-1] # Remove extra rows
    # y-axis
    fyvel=dt*(vel_f*yvel/dy)*(fpxy[:,2:]-fpxy[:,1:-1])
    fyvel=fyvel[1:-1,:] # Remove extra rows
    
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
    return fp

# pde function for no diffusion barrier
def pdef(wp,dt,dx,dy,xvel,yvel,vel_f,df,fp):

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
def pde(up,vp,wp,du,dv,dw,vel_u,vel_v,vel_w,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d):
    # Pass each array into their respective function
    up= pdef_db(wp,dt,dx,dy,xvel,yvel,vel_u,du,up)  
    vp= pdef_db(wp,dt,dx,dy,xvel,yvel,vel_v,dv,vp)
    wp= pdef(wp,dt,dx,dy,xvel,yvel,vel_w,dw,wp)
    
    # Apply the PDE to the new Arrays
    up = up + dt*(-a*up*vp)
    vp = vp + dt*(-b*up*vp+c*up*vp)
    wp = wp + dt*(d*up*vp)
    
    # Constant Conditions
    # Constant Brine Source in Center
    spotwidth=8 # This is half width in steps
    spotheight =1
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    vp[spotleft:spotright,ny-spotheight:ny]= 1      # Constant Brine Value
        
    # Return the resulting array
    return [up, vp, wp]