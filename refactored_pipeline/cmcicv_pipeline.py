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

import os
import pydicom
import matplotlib.pyplot as plt
import numpy as np




"""
==========================================================================
3D Image Reader

Read a series of images

Brendan Bouchard
20260126
Last Updated: 20260127
==========================================================================
"""

caseid = "001" # CHANGE THIS TO THE CURRENT CASE
fpath = "D:/cmcicv/mri_images/dicom/" # folder path of MRI images
pngoutpath = f"D:/cmcicv/mri_images/png/{caseid}/" # output path for .pngs
pixeloutpath = f"D:/cmcicv/mri_images/pixel_data/{caseid}/"
first_frame = 48 # first frame in sequence
last_frame = 52 # last frame in sequence

with os.scandir(fpath) as files:
    for i in files: # sequence through dicom files
        
        image_path = i.name 
        img_nopre = str.removeprefix(i.name, "IM_")
        step = str.removesuffix(img_nopre, ".dcm") # strips file name string to only the time step
        
        mri_image = os.path.join(fpath,image_path) # full path of dicom file
        
        # dicom file reader
        dicom_file = pydicom.dcmread(mri_image) # reads dicom file
        trigger_time = dicom_file.TriggerTime # Gets trigger time from metadata
        print("File ID: IM_00" + i.name + "\n" + "Trigger time: " + str(trigger_time) + "\n") # prints file ID and trigger time
       
        
        # plot result as image
        img = dicom_file.pixel_array # defines image from metadata
        plt.imshow(img, cmap=plt.cm.bone) # plots image on graph
        plt.title("Case ID: " + str(caseid) + " Timestep: " + step) 
        plt.savefig(pngoutpath + i.name + ".png")
        plt.show() # shows graph
        
        # saves pixel data to file
        pixels = dicom_file.pixel_array
        np.save(pixeloutpath + "IM_" + step + ".npy", pixels)
    

"""
==========================================================================
Volume Interpolator

Interpolates between points in the 

Brendan Bouchard
20260126
Last Updated: 20260126
==========================================================================
"""

