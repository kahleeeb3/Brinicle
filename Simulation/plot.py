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

#%% Functions
def func(x, a):
    return a * np.sqrt(x)
        
# (folder, tEnd, data2, L, dt, dtWindow, telap)
def plot(folder,data,param,time_vals):
    
    # Parse Input Data
    [Lx,Ly,res] = [param["Lx"], param["Ly"], param["res"]]
    [ tEnd , dt , dtWindow, time] = time_vals
    
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
    xdata = data[3][0,:t]
    ydata = data[3][1,:t]
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
    