# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import GMesh as gm
import os

cwd = os.getcwd()

####### DEV

# Escrito por Luís Cunha
# Última atualização : 12/02/2018


# Definição da malha
arquivo = "degrau"

print "Reading .msh file"
malha = gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN

# Parametros

nodes = len(x)
elenum = len(ien)
dt = 0.05
tempo = 150
Reynolds = 10
boundarynodes = len(malha.dirichlet_points)

# ---------------------------------------
# Wz, Psi e velocidade inicial
# ---------------------------------------

Psi = sp.zeros((nodes, tempo), dtype="float64")
Wz = sp.zeros((nodes, tempo), dtype="float64")
vx = sp.zeros(nodes, dtype="float64")
vy = sp.zeros(nodes, dtype="float64")


# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    K_local = sp.zeros((3, 3), dtype="float64")
    M_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    Gx_local = sp.zeros((3, 3), dtype="float64")
    Gy_local = sp.zeros((3, 3), dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    K_global = sp.zeros((_numnode, _numnode), dtype="float64")
    M_global = sp.zeros((_numnode, _numnode), dtype="float64")
    Gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    Gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

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
                K_local[i, j] = (b[i] * b[j] + c[i] * c[j]) / (4 * Area)
                Gx_local[i,j] = b[j] * (1/6.)
                Gy_local[i,j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                K_global[i_global, j_global] += K_local[i_local, j_local]
                M_global[i_global, j_global] += M_local[i_local, j_local]* (Area/12.)
                Gx_global[i_global, j_global] += Gx_local[i_local, j_local]
                Gy_global[i_global, j_global] += Gy_local[i_local, j_local]


    return  K_global, M_global, Gx_global, Gy_global

# POSSIVEL ERRO EM APPLY_DIRICHLET COM O ARRAY .DIRICHLET_POINTS CONSIDERANDO
# O MESMO PONTO MAIS DE UMA VEZ
def apply_dirichlet(bcpoints, nodes, _K):
    quantcc = len(bcpoints)
    KK = sp.copy(_K)
    BC = sp.zeros(nodes)

    for i in range(quantcc):
        index = int(bcpoints[i][0] - 1)
        value = bcpoints[i][1]
        for j in range(nodes):
            BC[j] -= value * _K[j, index]
            if j != index:
                KK[index, j] = 0
                KK[j, index] = 0
            else:
                KK[index, j] = 1
    return KK, BC

print "Assembling Matrices"

K, M, Gx, Gy = fem_matrix(x, y, elenum, nodes, ien)

Minv = sp.linalg.inv(M)
MinvLump = sp.zeros((nodes,nodes))

for i in range(nodes):
    for j in range(nodes):
        MinvLump[i, i] += Minv[i, j]

# ---------------------------------------
# Condições de contorno em K_psi
# ---------------------------------------

K_psi, ccpsi = apply_dirichlet(malha.dirichlet_points, nodes, K)

for i in range(boundarynodes):
    index = int(malha.dirichlet_points[i][0]-1)
    vy[index] = 0.0
    if malha.X[index] == 0.0:
        vx[index] = 1.0
    else:
        vx[index] = 0.0

Wz[0:nodes, 0] = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))

F_psi = sp.dot(M, Wz[0:nodes, 0]) + ccpsi
for i in range(boundarynodes):
    index = int(malha.dirichlet_points[i][0]-1)
    F_psi[index] = malha.dirichlet_points[i][1]

Psi[0:nodes, 0] = sp.linalg.solve(K_psi, F_psi)

#-------------------------------------------------
# Inicio de loop no tempo
#-------------------------------------------------

for t in range(0, tempo-1):
    BC_omega = sp.zeros((boundarynodes, 2), dtype="float64")

    print "Solving System " + str((float(t)/(tempo-1))*100) + "%"

    # B.C. Vorticidade

    Wcc = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))

    for i in range(boundarynodes):
        index = int(malha.dirichlet_points[i, 0])-1
        BC_omega[i][0] = index + 1
        BC_omega[i][1] = Wcc[index]
        vy[index] = 0.0
        if malha.X[index] == 0.0:
            vx[index] = 1.0
            BC_omega[i][1] = 0.0
        else:
            vx[index] = 0.0

    # Solução de Wz e Psi

    Conv = sp.dot(sp.diag(vx), Gx) + sp.dot(sp.diag(vy), Gy)
    LHS = (M/dt)+(K / Reynolds + Conv)
    LHS_omega, ccomega = apply_dirichlet(BC_omega, nodes, LHS)
    F_omega = sp.dot(M, Wz[0:nodes, t])/dt + ccomega / Reynolds

    for i in range(boundarynodes):
        index = int(malha.dirichlet_points[i, 0])-1
        F_omega[index] = Wcc[index]

    Wz[0:nodes, t+1] = sp.linalg.solve(LHS_omega, F_omega)

    F_psi = sp.dot(M, Wz[0:nodes, t+1]) + ccpsi
    for i in range(boundarynodes):
        index = int(malha.dirichlet_points[i][0]-1)
        F_psi[index] = malha.dirichlet_points[i][1]

    Psi[0:nodes, t+1] = sp.linalg.solve(K_psi, F_psi)

    vtk = io.InOut(x, y, ien, len(x), len(ien), Wz[:, t], Psi[:, t], vx, vy)
    vtk.saveVTK(cwd+"/results", arquivo + "-200-" + str(t+1))

    # Calculo de Vx e Vy

    vx = sp.dot(MinvLump, sp.dot(Gy, Psi[0:nodes, t+1]))
    vy = -1.0 * sp.dot(MinvLump, sp.dot(Gx, Psi[0:nodes, t+1]))

#-------------------------------------------------
# Fim de Loop no tempo
#-------------------------------------------------
