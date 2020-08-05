# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import os
from scipy import linalg
import InOut as Io

from assembly import fem_matrix as fem_py
from cyassembly import fem_matrix as fem_cy
import GMesh as gm
#import matplotlib.pyplot as plt
#from libmalha import lin2d, ret_tri


# L = 10.0
# N = 2
# x,y,ien = lin2d(N,L,N,L)
# ien = ret_tri(ien)
# nodes = len(x)
# elem = len(ien)

cwd = os.getcwd()

mesh_file = "viv.msh"
mesh = gm.GMesh(cwd + '/gTest/' + mesh_file)
x = mesh.X
y = mesh.Y
ien = mesh.IEN
nodes = len(x)
elem = len(ien)

K,M,Gx,Gy = fem_py(x, y, elem, nodes, ien)

#k2,m2,gx2,gy2 = fem_py(x, y, elem, nodes, ien)


psi = sp.zeros(nodes)

for i in range(nodes):
    psi[i] = y[i]

print('Solving Velocity')
#Minv = sp.linalg.inv(M)
# vx = sp.dot(Minv, sp.dot(Gy, Psi_new))
# vy = -1.0 * sp.dot(Minv, sp.dot(Gx, Psi_new))
print('vx')
vx = sp.linalg.solve(M, sp.dot(Gy, psi))
print('vy')
vy = -1.0 * sp.linalg.solve(M, sp.dot(Gx, psi))

print(vx[4244], vx[3415])
print(vy[4244], vy[3415])

print('Saving VTK')
vtk = Io.InOut(x, y, ien, len(x), len(ien), psi, None, None,
                None, None, vx, vy, None, None)
vtk.saveVTK(cwd+'/gTest', mesh_file)

exit()


# a = m
# a2 = m2
#
# print(a)
# print(type(a), type(a[-1,0]),'\n')
#
# print(a2)
# print(type(a2), type(a2[-1,0]),'\n')
#
# print(a2-a)

np.set_printoptions(suppress=True, precision = 4)

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow(k, cmap=plt.cm.Blues)
ax2.imshow(k2, cmap=plt.cm.Blues)

plt.figure(2)
plt.imshow(k2-k)
plt.colorbar()
plt.show()
