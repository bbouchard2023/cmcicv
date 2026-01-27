# -*- coding: utf-8 -*-
"""
==========================================================================
3D Image Reader

Read a series of images

Brendan Bouchard
20260126
Last Updated: 20260126
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
first_frame = 48 # first frame in sequence
last_frame = 52 # last frame in sequence

for i in range(first_frame,last_frame + 1): # sequence through dicom files
    
    image_path = "IM_00" + str(i) + ".dcm" 
    
    mri_image = os.path.join(fpath,image_path) # full path of dicom file
    
    # dicom file reader
    dicom_file = pydicom.dcmread(mri_image) # reads dicom file
    trigger_time = dicom_file.TriggerTime # Gets trigger time from metadata
    print("File ID: IM_00" + str(i) + "\n" + "Trigger time: " + str(trigger_time) + "\n") # prints file ID and trigger time
   
    
    # plot result as image
    img = dicom_file.pixel_array # defines image from metadata
    plt.imshow(img, cmap=plt.cm.bone) # plots image on graph
    plt.title("Case ID: " + str(caseid) + " Timestep: 00" + str(i)) 
    plt.savefig(pngoutpath + "IM_00" + str(i) + ".png")
    plt.show() # shows graph
    
    # saves pixel data to file
    pixels = dicom_file.pixel_array
    np.save(pixeloutpath + "IM_00" + str(i) + ".npy", pixels)
    
    