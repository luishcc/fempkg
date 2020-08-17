# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import GMesh as gm
import os

cwd = os.getcwd()

####### DEV

# Escrito por Luís Cunha
# Última atualização : 09/01/2018

# Código para a solução da EDP d²T               dT
#                              --- * alfa + Q = --- + v . grad(T)
#                              dx²               dt
# através de elementos Finitos.


# Definição da malha
arquivo = "teste"



print "Reading .msh file"
malha = gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN

# Parametros
nodes = len(x)
elenum = len(ien)
dt = 0.01
tempo = 600
teta = 1        # [0, 1] --> explicito - implicito


Alfa = sp.zeros(nodes) + 1.0

for i in range(nodes):
    if y[i] > 0.25:
        Alfa[i] = 0.1

# Montagem de matrizes

print "Assembling Matrices"
def fem_matrix(_x, _y, _numele, _numnode, _ien, _alfa):
    k_local = sp.zeros((3, 3), dtype="float64")
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    #gx_local = sp.zeros((3, 3), dtype="float64")
    #gy_local = sp.zeros((3, 3), dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    #gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    #gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):
        alfa = 0
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]
            alfa += _alfa[_ien[elem, i]]

        alfa = alfa * (1/3.0)
        a[0] = xx[0] * yy[1] - xx[1] * yy[0]
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        Area = (a[0] + a[1] + a[2]) / 2.

        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]

        for i in range(3):
            for j in range(3):
                k_local[i, j] = (b[i] * b[j] + c[i] * c[j]) / (4 * Area)
                #gx_local[i,j] = b[j] * (1/6.)
                #gy_local[i,j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]*alfa
                m_global[i_global, j_global] += m_local[i_local, j_local]* (Area/12.)
                #gx_global[i_global, j_global] += gx_local[i_local, j_local]
                #gy_global[i_global, j_global] += gy_local[i_local, j_local]


    return  k_global, m_global #, gx_global, gy_global


#K, M, Gx, Gy = fem_matrix(x, y, elenum, nodes, ien)

K, M = fem_matrix(x, y, elenum, nodes, ien, Alfa)


    # Implementando condições de contorno de DIRICHLET e INICIAL

print "Defining Boundary and Initial conditions"

quantcc = len(malha.dirichlet_points)
cc = sp.zeros(nodes)

len_neu = len(malha.Boundary_Neumann)
for i in range(len_neu):
    node1 = malha.Boundary_Neumann[i][0]+1
    node2 = malha.Boundary_Neumann[i][1]+1
    physgrp = malha.neumann_element_physgrp[i]
    for j in range(len(malha.neumann_points)):
        if node1 == int(malha.neumann_points[j, 0]) and physgrp == int(malha.neumann_points[j, 2]):
            value1 = malha.neumann_points[j, 1]
        if node2 == int(malha.neumann_points[j, 0]) and physgrp == int(malha.neumann_points[j, 2]):
            value2 = malha.neumann_points[j, 1]
    value = 0.5 * (value1 + value2)
    length = sp.sqrt((x[node1-1]-x[node2-1])**2 + (y[node1-1]-y[node2-1])**2)
    cc[node1-1] -= value*length
    cc[node2-1] -= value*length



#MQ = sp.dot(M, Q)
Mdt = M/dt
T_old = sp.zeros(nodes) + 1
LHS = Mdt + teta*(K)
KK = sp.copy(LHS)

for i in range(quantcc):
    index = int(malha.dirichlet_points[i][0]-1)
    value = malha.dirichlet_points[i][1]
    T_old[index] = malha.dirichlet_points[i, 1]
    for j in range(nodes):
        cc[j] -= value * LHS[j, index]
        if j != index:
            KK[index, j] = 0
            KK[j, index] = 0
        else:
            KK[index, j] = 1


    # Solucao do sistema

T_new = sp.zeros(nodes)

Matriz = (Mdt - (1 - teta)*(K))
for i in range(0, tempo):
    print "Solving System -- " + str(i+1) + " / "+str(tempo)

    B = sp.dot(Matriz, T_old) + cc
    for j in range(quantcc):
        index = int(malha.dirichlet_points[j][0]) - 1
        B[index] = malha.dirichlet_points[j][1]

    T_new = linalg.solve(KK, B)

    vtk = io.InOut(x, y, ien, len(x), len(ien), T_old, Alfa, None)
    vtk.saveVTK(cwd + "/results", arquivo + str(i))

    T_old = sp.copy(T_new)
