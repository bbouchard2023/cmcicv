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

def imgder(im1, im2, im3, delta):
    
    dx = np.zeros(2,2,2)
    dx[:,:,1] = np.array([-1, 1], [-1, 1])
    dx[:,:,2] = np.array([-1, 1], [-1, 1])
    dx = 0.25 * dx
    
    dy = np.zeros(2,2,2)
    dy[:,:,1] = np.array([-1 -1], [1, 1])
    dy[:,:,1] = np.array([-1, -1], [1, 1])
    dy = 0.25 * dy
    
    dz = np.zeros(2,2,2)
    dz[:,:,1] = np.array([-1, -1], [-1, -1])
    dz[:,:,2] = np.ones(2,2)
    dz = 0.25 * dz
    
    dt = np.ones(2,2,2)
    dt = 0.25 * dt
    
    Ix, Iy, Iz = np.gradient(im2, 1, 1, ((8 * 13) / (55 * 1.25))) 
    It = im3 - im1
    It = It / (1 * (2 * delta))
    
    return Ix, Iy, Iz, It

def opflow_lk(im1, im2, im3, radius, depth, delta):
    
    h1, w1, d1 = im1.size
    
    ux = np.zeros(im1.size)
    uy = ux
    uz = ux
    
    Ix, Iy, Iz, It = imgder(im1, im2, im3, delta)
    
    w = np.ones(len((np.arange(-radius,radius))),len((np.arange(-radius,radius))),len((np.arange(-radius,radius))))
    ww = w * w
    
    for i in range(radius + 1, h1 - radius):
        for j in range(radius + 1, w1 - radius):
            for k in range(depth + 1, d1 - depth):
                
                Ixsec = Ix[i - radius:i + radius, j - radius:j + radius, k - depth:k + depth]
                Iysec = Iy[i - radius:i + radius, j - radius:j + radius, k - depth:k + depth]
                Izsec = Iz[i - radius:i + radius, j - radius:j + radius, k - depth:k + depth]
                Itsec = It[i - radius:i + radius, j - radius:j + radius, k - depth:k + depth]
                
                A = np.zeros(3,3)
                B = np.zeros(3,1)
                
                A[1, 1] = sum(sum(sum(Ixsec^2 * ww)))
                A[1, 2] = sum(sum(sum(Ixsec * Iysec * ww)))
                A[1, 3] = sum(sum(sum(Ixsec * Izsec * ww)))
                
                A[2, 1] = sum(sum(sum(Iysec * Ixsec * ww)))
                A[2, 2] = sum(sum(sum(Iysec^2 * ww)))
                A[2, 3] = sum(sum(sum(Iysec * Izsec * ww)))
                
                A[3, 1] = sum(sum(sum(Izsec * Ixsec * ww)))
                A[3, 2] = sum(sum(sum(Izsec * Iysec * ww)))
                A[3, 3] = sum(sum(sum(Izsec^2 * ww)))
                
                B[1, 1] = sum(sum(sum(Ixsec * Itsec * ww)))
                B[2, 1] = sum(sum(sum(Iysec * Itsec * ww)))
                B[3, 1] = sum(sum(sum(Izsec * Itsec * ww)))
                
                pinvA = np.linalg.pinv(A)
                
                V = pinvA * -B # Calculates the flow vector
                ux[i, j, k] = V[1, 1]
                uy[i, j, k] = V[2, 1]
                uz[i, j, k] = V[3, 1]
                
                return ux, uy, uz

case = "001"
slice = 7
slicenum = str(slice)
mripath = "/home/brendan/cmcicv/mri_images/"
pxpath = mripath + f"pixel_data/{case}/"
pixelpath = natsorted(os.scandir(pxpath + f"slice_{slicenum}/"), key=lambda img: img.name)
radius = 10
depth = 5
delta = 0.0135

temparr = []
for pix in pixelpath:
    fname = pix.name
    arr = np.load(pxpath + f"slice_{slicenum}/" + fname)
    temparr.append(arr)

fullarr = np.stack(temparr,axis=-1)

for i in range(2, fullarr.shape[2]): 
    ux, uy, uz = opflow_lk(fullarr[:,:,i - 1], fullarr[:, :, i], fullarr[:, : , i + 1], radius, depth, delta) 



