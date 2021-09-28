from PIL import Image, ImageOps
import numpy as np #remove np
import cv2
import os

# Imports for Output Function
import math
import matplotlib.pyplot as plt
from matplotlib import cm

#%% Working
def image(name,x, y, angle):
    """
    Returns an array from an image file named image.png
    """
    
    im = Image.open(name)    #Open the image
    
    im = im.resize((y, x))  #Change Image Size
    im = im.rotate(angle)   #Rotate the Image
    im = im.convert('L')    #Convert to Black & White
    im = ImageOps.invert(im)   #Invert the image
    
    
    #im.save('image2.png') #Save the new image
    
    #Return the image as an array
    data = np.asarray(im)
    return data

def export(image_folder):
    """
    Converts image folder to a video
    """
    
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
    
def length(lp,wp):
    """
    For each row in the "ice array" find a value in which there is ice
    This should in essence give us the length of the ice in terms of pixels
    
    The pixels are plotted 90 degrees CCW which is (80,40)
    """
    #loop through the indecies of the middle
    center = int(len(wp)/2)
    for index, value in enumerate(wp[center,:]):
        if value > 0.01:
            print(index)
            length = 38 - index
            break
                
    return length

def improvedlength(lp,wp):
    length = len(wp)
    for x in range(length):
        #for each row: find an "ice value"
        for index, value in enumerate(wp[x,:]):
            if (value > 0.1) and (index < length):
                length = index
    return (len(wp[0,:])-length-2)/2
     
def plot(folder,tEnd, np,Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow):
    files = []   # This is for the output files.
    #windowmax = 30
    #res = 0.5
    t = 0        # This sets the first time point for the calculation.
    telap = 0    # This sets the time elapsed for the simulation.
    
    #text = 'u: Salt, v: Brine, w: Ice'
    #textx,texty,size = 10,-17,11
    
    while telap <= tEnd:    
    
        # Make the plots
        x = np.arange(0,Lx,dx) # Create the x data for plotting
        y = np.arange(0,Ly,dy) # Create the y data for plotting
        Y, X = np.meshgrid(y,x) # Create the X,Y matrices for plotting
        Z1 = u[:,:,t] # Pull the Z data matrix for plotting
        Z2 = v[:,:,t] # Pull the Z data matrix for plotting
        Z3 = w[:,:,t] # Pull the Z data matrix for plotting
        #Z4 = (u[:,:,t] + v[:,:,t]) # Pull the Z data matrix for plotting
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2) # Create the figure with subplots
        
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
        ax4.plot(l[0,:t],np.sqrt(l[0,:t]))
        #plt.xlim(0, 50)
        #plt.ylim(0, 50)
        #fig.colorbar(cf4, ax=ax4) # Add the colorbar
        ax4.set_xlabel('x (s)') # Label the x axis
        ax4.set_ylabel('y (mm)') # Label the y axis
        utitle = 'Length of w \n at {0:0.1f}s'.format(telap)
        ax4.set_title(utitle) # Title the plot
        #ax4.set_aspect(5) # Make the aspect ratio equal
        
        # plt.subplots_adjust(hspace=0.75,left=-0.05)
        #plt.text(textx, texty, text, fontsize=size)
        plt.tight_layout()
        # plt.subplots_adjust(left=-0.3)
    
        #plt.show() # This shows the plots as the code is running
        fname = 'images/Brinicle_%06d.png' % t # Create the file name for each plot
        print('Saving frame', fname) # Print the status update
        fig.savefig(fname, dpi=300) # Save the image
        files.append(fname) # Update the filename
        plt.close(fig) # Close the image so it doesn't show while the code is running
        plt.clf() # Clear the figrue to save memory
        t=t+math.floor(dtWindow/dt)            # Increment the storage counter
        telap=t*dt         # Increment the simulation time elapsed
        
def singlePlot(folder,tEnd, np,Lx, Ly, dx, dy, u, v, w, l, dt, dtWindow):
    files = []   # This is for the output files.
    t = math.floor(tEnd/dt)
    telap = tEnd
    
    #Make the plots
    x = np.arange(0,Lx,dx) # Create the x data for plotting
    y = np.arange(0,Ly,dy) # Create the y data for plotting
    Y, X = np.meshgrid(y,x) # Create the X,Y matrices for plotting
    Z1 = u[:,:,t] # Pull the Z data matrix for plotting
    Z2 = v[:,:,t] # Pull the Z data matrix for plotting
    Z3 = w[:,:,t] # Pull the Z data matrix for plotting
    #Z4 = (u[:,:,t] + v[:,:,t]) # Pull the Z data matrix for plotting
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2) # Create the figure with subplots
    
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
    ax4.plot(l[0,:t],np.sqrt(l[0,:t]))
    #plt.xlim(0, 50)
    #plt.ylim(0, 50)
    #fig.colorbar(cf4, ax=ax4) # Add the colorbar
    ax4.set_xlabel('x (s)') # Label the x axis
    ax4.set_ylabel('y (mm)') # Label the y axis
    utitle = 'Length of w \n at {0:0.1f}s'.format(telap)
    ax4.set_title(utitle) # Title the plot
    #ax4.set_aspect(5) # Make the aspect ratio equal
    
    # plt.subplots_adjust(hspace=0.75,left=-0.05)
    #plt.text(textx, texty, text, fontsize=size)
    plt.tight_layout()
    # plt.subplots_adjust(left=-0.3)

    #plt.show() # This shows the plots as the code is running
    fname = 'images/Brinicle_%06d.png' % t # Create the file name for each plot
    print('Saving frame', fname) # Print the status update
    fig.savefig(fname, dpi=300) # Save the image
    files.append(fname) # Update the filename
    plt.close(fig) # Close the image so it doesn't show while the code is running
    plt.clf() # Clear the figrue to save memory
    t=t+math.floor(dtWindow/dt)            # Increment the storage counter
    telap=t*dt         # Increment the simulation time elapsed
#%% Old PDE Function

#  Original
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
    spotwidth=2 # This is half width in steps
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    vp[spotleft:spotright,ny-5:ny-2]= 1      # Constant Brine Value
        
    # Return the resulting array
    return [up, vp, wp]

# Simplyfied
def ard_pdef(up,vp,wp,dt,dx,dy,xvel,yvel,vel_f,df,fp):
# This numerically integrates using finite differences and the forward
# Euler method.
        
    # Set up zero flux bc's in x.
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))

    # Make boundaries zero flux in y.
    fpxy = np.vstack((fpy[1,:], fpy, fpy[-2,:]))
            
    # Perform finite differences. (Diffusion)
    # Calculate differences
    fxx=dt*(df/dx**2)*(fpxy[:,2:]+fpxy[:,0:-2]-2*fpxy[:,1:-1])
    fxx=fxx[1:-1,:] # Remove extra rows
    # Calculate differences
    fyy=dt*(df/dy**2)*(fpxy[2:,:]+fpxy[0:-2,:]-2*fpxy[1:-1,:])
    fyy=fyy[:,1:-1] # Remove extra columns
   
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15;  # Set ff=1 to turn the fudge-factor off.
    
    # Calculate differences
    fxy=dt*((df/(dx**2+dy**2))*(fpxy[2:,2:]+fpxy[0:-2,0:-2]-2*fpxy[1:-1,1:-1]))/ff
    fyx=dt*((df/(dx**2+dy**2))*(fpxy[0:-2,2:]+fpxy[2:,0:-2]-2*fpxy[1:-1,1:-1]))/ff
    
    # Caclulate the advection along the x axis.
    fxvel=dt*(vel_f*xvel/dx)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fxvel=fxvel[:,1:-1] # Remove extra rows
    
    # Caclulate the advection along the y axis.
    #fyvel=dt*(vel_f*yvel/dy)*(fpxy[:,0:-2]-fpxy[:,1:-1]) # Upward Advection
    fyvel=dt*(vel_f*yvel/dy)*(fpxy[:,2:]-fpxy[:,1:-1]) # Downward Advection
    fyvel=fyvel[1:-1,:] # Remove extra rows
            
    # Combine diffusion and advection
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
        
    return fp

# Simplyfied (Calls ard_pdef())
def ard_pde1(up,vp,wp,du,dv,dw,vel_u,vel_v,vel_w,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d):
    # Pass each array into their respective function
    up2= ard_pdef(up,vp,wp,dt,dx,dy,xvel,yvel,vel_u,du,up)  
    vp2= ard_pdef(up,vp,wp,dt,dx,dy,xvel,yvel,vel_v,dv,vp)
    wp2= ard_pdef(up,vp,wp,dt,dx,dy,xvel,yvel,vel_w,dw,wp)
    
    # Apply the PDE to the new Arrays
    up = up2 + dt*(-a*up2*vp2)
    vp = vp2 + dt*(-b*up2*vp2+c*up2*vp2)
    wp = wp2 + dt*(d*up2*vp2)
    
    # Constant Conditions
    # Constant Brine Source in Center
    spotwidth=2 # This is half width in steps
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    vp[spotleft:spotright,ny-5:ny-2]= 1      # Constant Brine Value
        
    # Return the resulting array
    return [up, vp, wp]

# Adding diffusion barrier
def ard_pde2(up,vp,wp,du,dv,dw,vel_u,vel_v,vel_w,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d):
    # This numerically integrates using finite differences and the forward
    # Euler method.
    # Rates of diffusion matrix
    dum = du/(1+wp)
    dumy = np.column_stack((dum[:,1],dum,dum[:,-2]))
    dumxy = np.vstack((dumy[1,:], dumy, dumy[-2,:]))
    
    dvm = dv/(1+wp)
    dvmy = np.column_stack((dvm[:,1],dvm,dvm[:,-2]))
    dvmxy = np.vstack((dvmy[1,:], dvmy, dvmy[-2,:]))
        
    # zero flux boundary conditions
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
    upy = np.column_stack((up[:,1],up,up[:,-2]))
    vpy = np.column_stack((vp[:,1],vp,vp[:,-2]))
    wpy = np.column_stack((wp[:,1],wp,wp[:,-2]))
    
    upxy = np.vstack((upy[1,:], upy, upy[-2,:]))
    vpxy = np.vstack((vpy[1,:], vpy, vpy[-2,:]))
    wpxy = np.vstack((wpy[1,:], wpy, wpy[-2,:]))
    
    # Diffusion
    uxx=(dt/dx**2)*(dumxy[:,2:]*upxy[:,2:]+dumxy[:,0:-2]*upxy[:,0:-2]-2*dumxy[:,1:-1]*upxy[:,1:-1])
    uxx=uxx[1:-1,:] # Remove extra rows
    uyy=(dt/dy**2)*(dumxy[2:,:]*upxy[2:,:]+dumxy[0:-2,:]*upxy[0:-2,:]-2*dumxy[1:-1,:]*upxy[1:-1,:])
    uyy=uyy[:,1:-1] # Remove extra columns

    vxx=(dt/dx**2)*(dvmxy[:,2:]*vpxy[:,2:]+dvmxy[:,0:-2]*vpxy[:,0:-2]-2*dvmxy[:,1:-1]*vpxy[:,1:-1])
    vxx=vxx[1:-1,:] # Remove extra rows    
    vyy=(dt/dy**2)*(dvmxy[2:,:]*vpxy[2:,:]+dvmxy[0:-2,:]*vpxy[0:-2,:]-2*dvmxy[1:-1,:]*vpxy[1:-1,:])
    vyy=vyy[:,1:-1] # Remove extra columns
    
    """
    wxx=dt*(dw/dx**2)*(wpxy[:,2:]+wpxy[:,0:-2]-2*wpxy[:,1:-1])
    wxx=wxx[1:-1,:] # Remove extra rows
    wyy=dt*(dw/dy**2)*(wpxy[2:,:]+wpxy[0:-2,:]-2*wpxy[1:-1,:])
    wyy=wyy[:,1:-1] # Remove extra columns
    """
    wxx = 0     # no diffusion of ice 
    wyy = 0     # no diffusion of ice
    
    # Fudge Factor
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15;  # Set ff=1 to turn the fudge-factor off.
    
    # Diffusion (diagonal terms)
    uxy=((dt/(dx**2+dy**2))*(dumxy[2:,2:]*upxy[2:,2:]+dumxy[0:-2,0:-2]*upxy[0:-2,0:-2]-2*dumxy[1:-1,1:-1]*upxy[1:-1,1:-1]))/ff
    uyx=((dt/(dx**2+dy**2))*(dumxy[0:-2,2:]*upxy[0:-2,2:]+dumxy[2:,0:-2]*upxy[2:,0:-2]-2*dumxy[1:-1,1:-1]*upxy[1:-1,1:-1]))/ff
    
    vxy=((dt/(dx**2+dy**2))*(dvmxy[2:,2:]*vpxy[2:,2:]+dvmxy[0:-2,0:-2]*vpxy[0:-2,0:-2]-2*dvmxy[1:-1,1:-1]*vpxy[1:-1,1:-1]))/ff
    vyx=((dt/(dx**2+dy**2))*(dvmxy[0:-2,2:]*vpxy[0:-2,2:]+dvmxy[2:,0:-2]*vpxy[2:,0:-2]-2*dvmxy[1:-1,1:-1]*vpxy[1:-1,1:-1]))/ff
    
    wxy=dt*((dw/(dx**2+dy**2))*(wpxy[2:,2:]+wpxy[0:-2,0:-2]-2*wpxy[1:-1,1:-1]))/ff
    wyx=dt*((dw/(dx**2+dy**2))*(wpxy[0:-2,2:]+wpxy[2:,0:-2]-2*wpxy[1:-1,1:-1]))/ff
    
    # Advection
    # x-axis
    uxvel=dt*(vel_u*xvel/dx)*(upxy[0:-2,:]-upxy[1:-1,:])
    uxvel=uxvel[:,1:-1] # Remove extra rows
    vxvel=dt*(vel_v*xvel/dx)*(vpxy[0:-2,:]-vpxy[1:-1,:])
    vxvel=vxvel[:,1:-1] # Remove extra rows
    wxvel=dt*(vel_w*xvel/dx)*(wpxy[0:-2,:]-wpxy[1:-1,:])
    wxvel=wxvel[:,1:-1] # Remove extra rows
    
    # y-axis
    uyvel=dt*(vel_u*yvel/dy)*(upxy[:,0:-2]-upxy[:,1:-1])
    uyvel=uyvel[1:-1,:] # Remove extra rows
    vyvel=dt*(vel_v*yvel/dy)*(vpxy[:,2:]-vpxy[:,1:-1])
    vyvel=vyvel[1:-1,:] # Remove extra rows
    wyvel=dt*(vel_w*yvel/dy)*(wpxy[:,0:-2]-wpxy[:,1:-1])
    wyvel=wyvel[1:-1,:] # Remove extra rows
            
    # Combine diffusion with advection and kinetics using the forward Euler algorithm.
    up=up+(uxx+uyy+uxy+uyx)+(uxvel+uyvel) + dt*(-a*up*vp)     #(-k1*up*vp)
    vp=vp+(vxx+vyy+vxy+vyx)+(vxvel+vyvel) + dt*(-b*up*vp+c*up*vp)     #(-k1*up*vp)
    wp=wp+(wxx+wyy+wxy+wyx)+(wxvel+wyvel) + dt*(d*up*vp)     #(+k1*up*vp)
    
    # Constant Conditions
    # Constant Brine Source in Center
    spotwidth=2 # This is half width in steps
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    vp[spotleft:spotright,ny-5:ny-2]= 1      # Constant Brine Value
    
        
    return [up, vp, wp]

#%% New PDE Functions

# pde function that alters code for a diffusion barrier
def pdef1(up,vp,wp,dt,dx,dy,xvel,yvel,vel_f,df,fp):
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
def pdef2(up,vp,wp,dt,dx,dy,xvel,yvel,vel_f,df,fp):
    # zero flux boundary conditions
    # These pad the array with the second and second to last items,
    # Thus the last item has the same element on both sides.
    fpy = np.column_stack((fp[:,1],fp,fp[:,-2]))
    fpxy = np.vstack((fpy[1,:], fpy, fpy[-2,:]))
    
    # Diffusion
    fxx=dt*(df/dx**2)*(fpxy[:,2:]+fpxy[:,0:-2]-2*fpxy[:,1:-1])
    fxx=fxx[1:-1,:] # Remove extra rows
    fyy=dt*(df/dy**2)*(fpxy[2:,:]+fpxy[0:-2,:]-2*fpxy[1:-1,:])
    fyy=fyy[:,1:-1] # Remove extra columns
    
    # Fudge Factor
    # The included fudge-factor ff rounds out the square pixels
    # during diffusion. This factor downplays diagonal diffusion.
    ff=15  # Set ff=1 to turn the fudge-factor off.
    
    # Diffusion (diagonal terms)
    fxy=dt*((df/(dx**2+dy**2))*(fpxy[2:,2:]+fpxy[0:-2,0:-2]-2*fpxy[1:-1,1:-1]))/ff
    fyx=dt*((df/(dx**2+dy**2))*(fpxy[0:-2,2:]+fpxy[2:,0:-2]-2*fpxy[1:-1,1:-1]))/ff
    
    # Advection
    fxvel=dt*(vel_f*xvel/dx)*(fpxy[0:-2,:]-fpxy[1:-1,:])
    fxvel=fxvel[:,1:-1] # Remove extra rows
    fyvel=dt*(vel_f*yvel/dy)*(fpxy[:,0:-2]-fpxy[:,1:-1])
    fyvel=fyvel[1:-1,:] # Remove extra rows
    
    fp=fp+(fxx+fyy+fxy+fyx)+(fxvel+fyvel)
    
    return fp

# Calls pdef1() and pdef2() above
def pde(up,vp,wp,du,dv,dw,vel_u,vel_v,vel_w,xvel,yvel,nx,ny,dt,dx,dy,a,b,c,d):
    # Pass each array into their respective function
    up2= pdef1(up,vp,wp,dt,dx,dy,xvel,yvel,vel_u,du,up)  
    vp2= pdef1(up,vp,wp,dt,dx,dy,xvel,yvel,vel_v,dv,vp)
    wp2= pdef2(up,vp,wp,dt,dx,dy,xvel,yvel,vel_w,dw,wp)
    
    # Apply the PDE to the new Arrays
    up = up2 + dt*(-a*up2*vp2)
    vp = vp2 + dt*(-b*up2*vp2+c*up2*vp2)
    wp = wp2 + dt*(d*up2*vp2)
    
    # Constant Conditions
    # Constant Brine Source in Center
    spotwidth=2 # This is half width in steps
    spotleft=int(np.round(nx/2))-spotwidth   # Determine the left edge
    spotright=int(np.round(nx/2))+spotwidth  # Determine the right edge
    vp[spotleft:spotright,ny-5:ny-2]= 1      # Constant Brine Value
        
    # Return the resulting array
    return [up, vp, wp]