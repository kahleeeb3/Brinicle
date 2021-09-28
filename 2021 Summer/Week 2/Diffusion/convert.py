from PIL import Image, ImageOps
import numpy
import cv2
import os

def image():
    im = Image.open("image.png")
    #Make the new image half the width and half the height of the original image
    print('Image Size:',im.size)
    im = im.rotate(180)
    print('Image Size:',im.size)
    im = im.resize((200, 150))
    im = im.convert('L')
    #im = ImageOps.invert(im)
    #Save the cropped image
    im.save('image2.png')
    
    data = numpy.asarray(im)
    print("test")
    return data

def export(image_folder):
    print("exporting to video")
    video_name = 'video.avi'
    fps = 10

    
    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape
    
    video = cv2.VideoWriter(video_name, 0, fps, (width,height))
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))
    
    cv2.destroyAllWindows()
    video.release()
    print("Completed")