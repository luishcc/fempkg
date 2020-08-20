# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import MESH as ms
from apply_bc import apply_bc_dirichlet as bc_ap
from cython_func.assembly import ca
from timeit import default_timer as timer



file = "viv3.msh"

mesh = ms.Mesh(file)

y = mesh.y
x = mesh.x
ien = mesh.ien

k,m,gx,gy = ca.fem_matrix(mesh)

LHS = k + m

#--------------------------
Boundary = np.zeros(mesh.num_nodes)
lista = []

bc = np.zeros(mesh.num_nodes)
print(len(bc))
cylinder = mesh.get_boundary_with_name('cylinder')

for i in range(mesh.num_nodes):
    if y[i] == 0.0 :
        Boundary[i] = i
        bc[i] = 0
    elif y[i] == 10.0 :
        Boundary[i] = i
        bc[i] = 10
    elif i in cylinder:
        Boundary[i] = i
        bc[i] = 5
    elif x[i] == 0.0:
        bc[i] = 10 * mesh.y[i]
        Boundary[i] = i
    elif x[i] == 32.5:
        Boundary[i] = i
    else:
        lista.append(i)

Boundary = np.delete(Boundary, lista, axis=0)
num_bc = len(Boundary)


t1 = timer()
l1 = np.copy(LHS)
ccpsi = np.zeros(mesh.num_nodes)

for i in range(num_bc):
    index = int(Boundary[i])
    value = bc[index]
    if y[index] == 0.0 or y[index] == 10.0 or \
    (index in cylinder) or x[index] == 0:
        for j in range(mesh.num_nodes):
            ccpsi[index] -= value * LHS[j, index]
            if j != index:
                l1[index, j] = 0
                l1[j, index] = 0
            else:
                l1[index, j] = 1

t2 = timer()
#-----------------------------------------------

#mesh.set_neighbours('vivC')
mesh.set_neighbours()
mesh.set_dirichlet_nodes(['inlet', 'top', 'bot', 'cylinder'])
print(mesh.dirichlet_nodes)
print(Boundary)

t3 = timer()
l2, cc2 = bc_ap(LHS, mesh, bc)
t4 = timer()

print(t2-t1)
print(t4-t3)


np.set_printoptions(suppress=True, precision = 4)

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow(l1, cmap=plt.cm.Blues)
ax2.imshow(l2, cmap=plt.cm.Blues)

plt.figure(2)
plt.imshow(l2-l1)
plt.colorbar()
plt.show()
