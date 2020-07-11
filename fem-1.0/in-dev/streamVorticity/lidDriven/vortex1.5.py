# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import os

cwd = os.getcwd()

arquivo = "vortice"

print "Reading .msh file"
malha = Gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN
nodes = len(x)
num_ele = len(ien)

dt = 0.05
tempo = 1000
Re = 10



# ---------------------------------------
# Wz, Psi e velocidade inicial
# ---------------------------------------

Psi_new = sp.zeros(nodes, dtype="float64")
Wz_new = sp.zeros(nodes, dtype="float64")
vx = sp.zeros(nodes, dtype="float64")
vy = sp.zeros(nodes, dtype="float64")

# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_local = sp.zeros((3, 3), dtype="float64")
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    gx_local = sp.zeros((3, 3), dtype="float64")
    gy_local = sp.zeros((3, 3), dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]

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
                gx_local[i,j] = b[j] * (1/6.)
                gy_local[i,j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]
                m_global[i_global, j_global] += m_local[i_local, j_local]* (Area/12.)
                gx_global[i_global, j_global] += gx_local[i_local, j_local]
                gy_global[i_global, j_global] += gy_local[i_local, j_local]


    return  k_global, m_global, gx_global, gy_global

print "Assembling Matrices"
K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

K_Re = K/Re
Mdt = M/dt
Mlump = sp.zeros((nodes,nodes))

Boundary = sp.zeros(nodes)
lista = []
for i in range(nodes):
    if x[i] == 0.0 or y[i] == 0.0 or y[i] == 1.0 or x[i] == 1.0:
        Boundary[i] = i
    else:
        lista.append(i)
    for j in range(nodes):
        Mlump[i, i] += M[i, j]



MinvLump = linalg.inv(Mlump)
Boundary = sp.delete(Boundary, lista, axis=0)
num_bc = len(Boundary)



# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------

print "Initial and Boundary Conditions"

for i in Boundary:
    j = int(i)
    vy[j] = 0.0
    if y[j] == 1.0:
        vx[j] = 1.0
    else:
        vx[j] = 0.0

Wz_old = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))

K_psi = sp.copy(K)
ccpsi = sp.zeros(nodes)

for i in range(num_bc):
    index = int(Boundary[i])
    value = 0.0
    for j in range(nodes):
        ccpsi[j] -= value * K[j, index]
        if j != index:
            K_psi[index, j] = 0
            K_psi[j, index] = 0
        else:
            K_psi[index, j] = 1

F_psi = sp.dot(M, Wz_old) + ccpsi
for i in range(num_bc):
    index = int(Boundary[i])
    F_psi[index] = 0.0

Psi_old = sp.linalg.solve(K_psi, F_psi)

# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

for t in range(0, tempo):
    print "Iteration: ", t, "/" + str(tempo-1)

    for i in Boundary:
        j = int(i)
        if y[j] == 1.0:
            vx[j] = 1.0
	if y[j] == 0.0:
	    vx[j] = 0.0       
	else:
            vy[j] = 0.0

    # B.C. Vorticidade
    Wcc = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))
    ccomega = sp.zeros(nodes)

    # Solução de Wz e Psi
    Conv = vx*Gx + vy*Gy
#    Conv = sp.dot(sp.diag(vx), Gx) + sp.dot(sp.diag(vx), Gx)
    LHS = Mdt + K_Re + Conv
    LHS_omega = sp.copy(LHS)

    for i in range(num_bc):
        index = int(Boundary[i])
        value = Wcc[index]
        for j in range(nodes):
            ccomega[j] -= value * LHS[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    F_omega = sp.dot(Wz_old, Mdt) + ccomega / float(Re) #- sp.dot(Conv, Wz_old)

    for i in range(num_bc):
        index = int(Boundary[i])
        F_omega[index] = Wcc[index]

    Wz_new = sp.linalg.solve(LHS_omega, F_omega)

    F_psi = sp.dot(M, Wz_new) + ccpsi
    for i in range(num_bc):
        index = int(Boundary[i])
        F_psi[index] = 0.0

    Psi_new = sp.linalg.solve(K_psi, F_psi)

    vtk = Io.InOut(x, y, ien, len(x), len(ien), Wz_old, Psi_old, vx, vy)
    vtk.saveVTK(cwd+"/results", arquivo + "Re"+ str(Re) +"-"+ str(t+1))

    Psi_old = sp.copy(Psi_new)
    Wz_old = sp.copy(Wz_new)

    # Calculo de Vx e Vy

    vx = sp.linalg.solve(M, sp.dot(Gy,Psi_new))
    vy = -1.0 * sp.linalg.solve(M, sp.dot(Gx,Psi_new))

#----------------- Fim de Loop -------------------
#-------------------------------------------------

print "Done - .vtk saved in /results"
