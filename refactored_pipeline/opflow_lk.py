#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
Optical Flow for MRI Chest Displacement

Hopefully will allow us to measure chest displacement to mm accuracy using Lucas-Kanade algorithm and optical flow

Brendan Bouchard
20260319
Last Updated: 20260319
==========================================================================
"""
import os
from natsort import natsorted

rootpath = "/home/brendan/cmcicv/mri_images/dicom/"

mripath = "/home/brendan/cmcicv/mri_images/"

try:
    os.mkdir("/home/brendan/cmcicv/mri_images/dicom_slices")
    outpath = "/home/brendan/cmcicv/mri_images/dicom_slices/"
except FileExistsError:
    print("Path already exists")
    outpath = "/home/brendan/cmcicv/mri_images/dicom_slices/"

sorteddir = natsorted(os.scandir(rootpath), key=lambda img: img.name)
slice_num = 1


for imgs in sorteddir:
    print(imgs.name)
    fname = imgs.name
    rpre = str.removeprefix(fname, "IM_")
    rsuf = str.removesuffix(rpre, ".dcm")
    num = str.lstrip(rsuf, "0")
    num = int(num)
    num = num - 48
    print(num)
    
    if num % 30 == 0:
        print("new slice")
        
    
    