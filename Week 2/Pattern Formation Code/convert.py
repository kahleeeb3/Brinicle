from PIL import Image, ImageOps
import numpy
import cv2
import os

# Imports for Output Function
import math
import matplotlib.pyplot as plt
from matplotlib import cm

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
    data = numpy.asarray(im)
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
    

def plot(folder,tEnd,displayUandV, np,Lx, Ly, dx, dy, u, v, dt, dtWindow):
    files = []   # This is for the output files.
    t = 0        # This sets the first time point for the calculation.
    telap = 0    # This sets the time elapsed for the simulation.
    plot_u_min = 0
    plot_u_max = 100
    plot_v_min = 0
    plot_v_max = 100
    
    while telap <= tEnd:    
    
        # Make the plots
        if displayUandV: # If the Boolean is set to plot both u and v
            x = np.arange(0,Lx,dx) # Create the x data for plotting
            y = np.arange(0,Ly,dy) # Create the y data for plotting
            X, Y = np.meshgrid(y,x) # Create the X,Y matrices for plotting
            Z1 = u[:,:,t] # Pull the Z data matrix for plotting
            Z2 = v[:,:,t] # Pull the Z data matrix for plotting
            fig, (ax1, ax2) = plt.subplots(2,1) # Create the figure with subplots
            # Create the filled countour plot with colormap and manual levels
            #cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm,levels=np.arange(plot_u_min,plot_u_max,0.05))
            cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm)
            fig.colorbar(cf1, ax=ax1) # Add the colorbar
            ax1.set_xlabel('x (mm)') # Label the x axis
            ax1.set_ylabel('y (mm)') # Label the y axis
            utitle = 'Concentration of u at {0:.1f}s'.format(telap)
            ax1.set_title(utitle) # Title the plot
            ax1.set_aspect('equal') # Make the aspect ratio equal
            # Create the filled countour plot with colormap and manual levels
            #cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm,levels=np.arange(plot_v_min,plot_v_max,0.05))
            cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm)
            fig.colorbar(cf2, ax=ax2) # Add the colorbar
            ax2.set_xlabel('x (mm)') # Label the x axis
            ax2.set_ylabel('y (mm)') # Label the y axis
            vtitle = 'Concentration of v at {0:.1f}s'.format(telap)
            ax2.set_title(vtitle) # Title the plot
            ax2.set_aspect('equal') # Make the aspect ratio equal
            plt.subplots_adjust(hspace=0.75,left=-0.5)
        fname = f'{folder}/Patterns_%06d.png' % t # Create the file name for each plot
        print('Saving frame', fname) # Print the status update
        fig.savefig(fname, dpi=300) # Save the image
        files.append(fname) # Update the filename
        plt.close(fig) # Close the image so it doesn't show while the code is running
        t=t+math.floor(dtWindow/dt)            # Increment the storage counter
        telap=t*dt         # Increment the simulation time elapsed   
        
def bruss(folder,tEnd,displayUandV, np,Lx, Ly, dx, dy, u, v, dt, dtWindow):
    files = []   # This is for the output files.
    t = 0        # This sets the first time point for the calculation.
    telap = 0    # This sets the time elapsed for the simulation.
    plot_u_min = 0
    plot_u_max = 6.4
    plot_v_min = 0
    plot_v_max = 6.4
    
    while telap <= tEnd:    
    
        # Make the plots
        if displayUandV: # If the Boolean is set to plot both u and v
            x = np.arange(0,Lx,dx) # Create the x data for plotting
            y = np.arange(0,Ly,dy) # Create the y data for plotting
            X, Y = np.meshgrid(y,x) # Create the X,Y matrices for plotting
            Z1 = u[:,:,t] # Pull the Z data matrix for plotting
            Z2 = v[:,:,t] # Pull the Z data matrix for plotting
            fig, (ax1, ax2) = plt.subplots(2,1) # Create the figure with subplots
            # Create the filled countour plot with colormap and manual levels
            cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm,levels=np.arange(plot_u_min,plot_u_max,0.05))
            # cf1 = ax1.contourf(X,Y,Z1,cmap=cm.coolwarm)
            fig.colorbar(cf1, ax=ax1) # Add the colorbar
            ax1.set_xlabel('x (mm)') # Label the x axis
            ax1.set_ylabel('y (mm)') # Label the y axis
            utitle = 'Concentration of u at {0:.1f}s'.format(telap)
            ax1.set_title(utitle) # Title the plot
            ax1.set_aspect('equal') # Make the aspect ratio equal
            # Create the filled countour plot with colormap and manual levels
            cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm,levels=np.arange(plot_v_min,plot_v_max,0.05))
            # cf2 = ax2.contourf(X,Y,Z2,cmap=cm.coolwarm)
            fig.colorbar(cf2, ax=ax2) # Add the colorbar
            ax2.set_xlabel('x (mm)') # Label the x axis
            ax2.set_ylabel('y (mm)') # Label the y axis
            vtitle = 'Concentration of v at {0:.1f}s'.format(telap)
            ax2.set_title(vtitle) # Title the plot
            ax2.set_aspect('equal') # Make the aspect ratio equal
            plt.subplots_adjust(hspace=0.75,left=-0.5)
        else: # If the Boolean is set to plot only v
            x = np.arange(0,Lx,dx) # Create the x data for plotting
            y = np.arange(0,Ly,dy) # Create the y data for plotting
            X, Y = np.meshgrid(y,x) # Create the X,Y matrices for plotting
            Z = v[:,:,t] # Pull the Z data matrix for plotting
            fig, ax = plt.subplots() # Create the figure with subplots
            # Create the filled countour plot with colormap and manual levels
            cf = ax.contourf(X,Y,Z,cmap=cm.coolwarm,levels=np.arange(plot_v_min,plot_v_max,0.05))
            fig.colorbar(cf, ax=ax) # Add the colorbar
            ax.set_xlabel('x (mm)') # Label the x axis
            ax.set_ylabel('y (mm)') # Label the y axis
            vtitle = 'Concentration of v at {0:.1f}s'.format(telap)
            ax.set_title(vtitle) # Title the plot
            ax.set_aspect('equal') # Make the aspect ratio equal 
        # plt.show() # This shows the plots as the code is running
        fname = f'{folder}/Patterns_%06d.png' % t # Create the file name for each plot
        print('Saving frame', fname) # Print the status update
        fig.savefig(fname, dpi=300) # Save the image
        files.append(fname) # Update the filename
        plt.close(fig) # Close the image so it doesn't show while the code is running
        t=t+math.floor(dtWindow/dt)            # Increment the storage counter
        telap=t*dt         # Increment the simulation time elapsed   