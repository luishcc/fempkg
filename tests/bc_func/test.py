# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import MESH as ms
import GMesh as gm
from cython_func.assembly import ca
from cython_func.assembly import assembly as a



file = "1valid.msh"

mesh1 = gm.GMesh(file)
mesh2 = ms.Mesh(file)

k,m,gx,gy = a.fem_matrix(mesh1.x, mesh1.y, len(mesh1.ien), mesh1.Nnodes,
                            mesh1.ien)
k2,m2,gx2,gy2 = ca.fem_matrix(mesh2)

np.set_printoptions(suppress=True, precision = 4)

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow(m, cmap=plt.cm.Blues)
ax2.imshow(m2, cmap=plt.cm.Blues)

plt.figure(2)
plt.imshow(m2-m)
plt.colorbar()
plt.show()
