#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
Optical Flow for MRI Chest Displacement

Hopefully will allow us to measure chest displacement to mm accuracy using Lucas-Kanade algorithm and optical flow

Brendan Bouchard
20260319
Last Updated: 20260401
==========================================================================
"""
import os
from natsort import natsorted
import time

case = "001"

rootpath = "/home/brendan/cmcicv/mri_images/dicom/"

mripath = "/home/brendan/cmcicv/mri_images/"

pxpath = mripath + f"pixel_data/{case}/"

try:
    os.mkdir("/home/brendan/cmcicv/mri_images/dicom_slices")
    outpath = "/home/brendan/cmcicv/mri_images/dicom_slices/"
except FileExistsError:
    print("Path already exists")
    outpath = "/home/brendan/cmcicv/mri_images/dicom_slices/"

sorteddir = natsorted(os.scandir(rootpath), key=lambda img: img.name)
slice_num = 1


for imgs in sorteddir:
    fname = imgs.name
    print(fname)
    rpre = str.removeprefix(fname, "IM_")
    rsuf = str.removesuffix(rpre, ".dcm")
    num = str.lstrip(rsuf, "0")
    num = int(num)
    num = num - 48
    if num % 30 == 0 and num != 0:
       print("new slice")
       slice_num += 1
    try:
        os.mkdir(outpath + f"slice_{slice_num}/")
        slicepath = outpath + f"slice_{slice_num}/"
    except FileExistsError:
        slicepath = outpath + f"slice_{slice_num}/"
    
    time.sleep(0.4)
    fpath = slicepath + fname
    print(fpath)
    
    with open(fpath, 'w') as f:
        f.write(f"Saving {imgs} to {slicepath}")

# %%
"""
==========================================================================
Optical Flow for MRI Chest Displacement

Hopefully will allow us to measure chest displacement to mm accuracy using Lucas-Kanade algorithm and optical flow

Brendan Bouchard
20260402
Last Updated: 20260402
==========================================================================
"""

import os
from natsort import natsorted
import numpy as np

temparr = []
slice = 7

slicenum = str(slice)

pixelpath = natsorted(os.scandir(pxpath + f"slice_{slicenum}/"), key=lambda img: img.name)

for pix in pixelpath:
    fname = pix.name
    arr = np.load(pxpath + f"slice_{slicenum}/" + fname)
    temparr.append(arr)

fullarr = np.stack(temparr,axis=-1)

for i in range(1, len(fullarr))