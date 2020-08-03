# -*- coding: utf-8 -*-
import numpy as np
from assembly import fem_matrix as fem_py
from cyassembly import fem_matrix as fem_cy
import GMesh as gm
import matplotlib.pyplot as plt
from libmalha import lin2d, ret_tri


# L = 10.0
# N = 2
# x,y,ien = lin2d(N,L,N,L)
# ien = ret_tri(ien)
# nodes = len(x)
# elem = len(ien)

mesh_file = "1valid.msh"
mesh = gm.GMesh(mesh_file)
x = mesh.X
y = mesh.Y
ien = mesh.IEN
nodes = len(x)
elem = len(ien)

k,m,gx,gy = fem_py(x, y, elem, nodes, ien)

k2,m2,gx2,gy2 = fem_cy(x, y, elem, nodes, ien)

np.set_printoptions(suppress=True, precision = 4)

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

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow(k, cmap=plt.cm.Blues)
ax2.imshow(k2, cmap=plt.cm.Blues)

plt.figure(2)
plt.imshow(k2-k)
plt.colorbar()
plt.show()
