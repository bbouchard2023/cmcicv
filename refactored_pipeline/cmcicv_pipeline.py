# -*- coding: utf-8 -*-
"""
==========================================================================
CMCICV Pipeline

Takes MRI images and creates a geometry for displacement evaluation

Brendan Bouchard
20260126
Last Updated: 20260127
==========================================================================
"""



# %%


"""
==========================================================================
3D Image Reader

Read a series of images

Brendan Bouchard
20260126
Last Updated: 20260128
==========================================================================
"""

import os
import pydicom
import matplotlib.pyplot as plt
import numpy as np



caseid = "001" # CHANGE THIS TO THE CURRENT CASE
fpath = "D:/cmcicv/mri_images/dicom/" # folder path of MRI images
pngoutpath = f"D:/cmcicv/mri_images/png/{caseid}/" # output path for .pngs
pixeloutpath = f"D:/cmcicv/mri_images/pixel_data/{caseid}/"
img_slice = 0

with os.scandir(fpath) as files:
    for i in files: # sequence through dicom files
        
        image_path = i.name 
        img_nopre = str.removeprefix(i.name, "IM_")
        step = str.removesuffix(img_nopre, ".dcm") # strips file name string to only the time step
        
        mri_image = os.path.join(fpath,image_path) # full path of dicom file
        
        # dicom file reader
        dicom_file = pydicom.dcmread(mri_image) # reads dicom file
        trigger_time = dicom_file.TriggerTime # Gets trigger time from metadata
        print("File ID: IM_" + step + "\n" + "Trigger time: " + str(trigger_time) + "\n") # prints file ID and trigger time
       
        if (int(step) - 48) % 30 == 0:
           img_slice += 1
           os.mkdir(pngoutpath + f"slice_{img_slice}") 
           slice_path = pngoutpath + f"slice_{img_slice}/"
       
        # plot result as image
        img = dicom_file.pixel_array # defines image from metadata
        plt.imshow(img, cmap=plt.cm.bone) # plots image on graph
        plt.title("Case ID: " + str(caseid) + " Timestep: " + step) 
        plt.savefig(slice_path + "IM_" + step + ".png")
        plt.show() # shows graph
        
        
        
        # saves pixel data to file
        pixels = dicom_file.pixel_array
        np.save(pixeloutpath + "IM_" + step + ".npy", pixels)
    
# %%


"""
==========================================================================
Volume Structuring

Restructures data to temporal volume of slices

Brendan Bouchard
20260128
Last Updated: 20260128
==========================================================================
"""

import os
import pydicom
import matplotlib.pyplot as plt
import numpy as np



caseid = "001" # CHANGE THIS TO THE CURRENT CASE
pixeloutpath = f"D:/cmcicv/mri_images/pixel_data/{caseid}/"
voloutpath = f"D:/cmcicv/mri_images/volume_data/{caseid}/"
img_slice = 0
frames = []


with os.scandir(pixeloutpath) as files:
    for i in files:
        
        data = np.load(i.path) # imports data from each .npy file
        img = i.name # gets file name
        img_nopre = str.removeprefix(i.name, "IM_") 
        img_nosfx = str.removesuffix(img_nopre, ".npy")
        step = int(img_nosfx) # strips file name down to step number
 
        frames.append(data) # creates a list of all the arrays

        

        
        if (step - 48) % 30 == 0:
            img_slice += 1 # increments for each slice

    
stored_data = np.stack(frames,axis=2) # stores 2D data on top of each other (along the third dimension)       

mri_slices = np.array_split(stored_data,img_slice,axis=2) # splits each slice (30 frames per slice)

for idx, mrislice in enumerate(mri_slices, start=1):
    np.save(voloutpath + f"slice_{idx}.npy", mrislice)


# %%

"""
==========================================================================
Video Creator

Takes images from each slice and creates a video; image size is 256x256

Brendan Bouchard
20260128
Last Updated: 20260128
==========================================================================
"""

import os
import pydicom
import matplotlib.pyplot as plt
import numpy as np
import cv2 as cv



caseid = "001" # CHANGE THIS TO THE CURRENT CASE
pngoutpath = f"D:/cmcicv/mri_images/png/{caseid}/" # output path for .pngs
vidoutpath = f"D:/cmcicv/mri_images/videos/{caseid}/" # output path for .mp4s
fps = 1.5
imgslice = 0
numdir = len([name for name in os.listdir(pngoutpath)
              if os.path.isdir(pngoutpath + name)])
fps = 1.5


for i in range(1,numdir + 1):

# NOT WORKING
    vid_name = vidoutpath + f"slice_{i}.avi"
    fourcc = cv.VideoWriter_fourcc(*"MJPG")
    vid = cv.VideoWriter(vid_name, fourcc, fps, (256, 256))
    with os.scandir(pngoutpath + f"slice_{i}") as files:
        for img in files:
            
            images = cv.imread(img.path, cv.IMREAD_GRAYSCALE)
            images = cv.normalize(images, None, 0, 255, cv.NORM_MINMAX).astype(np.uint8)
            images = cv.cvtColor(images, cv.COLOR_GRAY2BGR)
            
            vid.write(images)
            
        vid.release()
cv.destroyAllWindows()
        
