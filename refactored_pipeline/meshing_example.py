# -*- coding: utf-8 -*-
"""
==========================================================================
Meshing Example

Creates a box then exports it as a .msh

Brendan Bouchard
20260309
Last Updated: 20260309
==========================================================================
"""

import numpy as np
import meshio

l0 = 1.0
A = 0.2
omega = 2*np.pi
t = 0.1
frames = 20

def box_size(t):
    return l0 + (A * np.sin(omega * t))

x = np.linspace(0,1,10)
y = np.linspace(0,1,10)
z = np.linspace(0,1,10)

u,v,w = np.meshgrid(x,y,z,indexing='ij')
nodes = np.vstack([u.ravel(),v.ravel(),w.ravel()]).T

def deform_nodes(nodes, scale):
    center = np.array([0.5,0.5,0.5])
    shifted = nodes - center
    scaled = shifted * scale
    return scaled + center


scale = box_size(t)


hex_cells = []

lenx = len(x)
leny = len(y)
lenz = len(z)

for i in range(lenx-1):
    for j in range(leny-1):
        for k in range(lenz-1):
            n0 = i*leny*lenz + j*lenz + k
            n1 = (i+1)*leny*lenz + j*lenz + k
            n2 = (i+1)*leny*lenz + (j+1)*lenz + k
            n3 = i*leny*lenz + (j+1)*lenz + k
            n4 = i*leny*lenz + j*lenz + (k+1)
            n5 = (i+1)*leny*lenz + j*lenz + (k+1)
            n6 = (i+1)*leny*lenz + (j+1)*lenz + (k+1)
            n7 = i*leny*lenz + (j+1)*lenz + (k+1)

            hex_cells.append([n0,n1,n2,n3,n4,n5,n6,n7])
hex_cells = np.array(hex_cells) + 1

scales = np.linspace(0.8,1.2,frames//2)
scales = np.concatenate([scales, scales[::-1]])

for t, scale in enumerate(scales):
    new_nodes = deform_nodes(nodes, scale)
    msh = meshio.Mesh(
        points=new_nodes,
        cells=[("hexahedron", hex_cells)]
        )
    fid = f"D:/cmcicv/sim_data_processing/mesh_creation/box_{t:02d}.msh"
    msh.write(fid, file_format="gmsh22")