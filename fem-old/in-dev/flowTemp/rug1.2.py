# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import os

cwd = os.getcwd()

#pra rugosidade mudar x - 5.5

arquivo = "rugosidade"

malha = Gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN
nodes = len(x)
num_ele = len(ien)

dt = 0.1
tempo = 120
Re = 10.0
Pr = 1.0

# ---------------------------------------
# Wz, Psi e velocidade inicial
# ---------------------------------------

Psi_new = sp.zeros(nodes, dtype="float64")
Wz_new = sp.zeros(nodes, dtype="float64")
Temp_old = sp.zeros(nodes, dtype="float64")
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

K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

K_Re = K/Re
Mdt = M/dt
MLump = sp.zeros((nodes,nodes))

lista = []
velx_nodes= sp.zeros(nodes)
for i in range(nodes):
    if x[i] == 0 and (1.0> y[i] >0.0):
        vx[i] = 1.0
        velx_nodes[i] = i
    else:
        lista.append(i)
    for j in range(nodes):
        MLump[i, i] += M[i, j]

velx_nodes = sp.delete(velx_nodes, lista, axis=0)
velx_quant = len(velx_nodes)

MinvLump = linalg.inv(MLump)


# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------
psi_dirichlet = sp.zeros(nodes)
bc_omega = sp.zeros(nodes)

psi_boundary = len(malha.dirichlet_points)
for i in range(psi_boundary):
    index = int(malha.dirichlet_points[i][0]-1)
    value = malha.dirichlet_points[i][1]
    psi_dirichlet[index] = value

Wz_old = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))

K_psi = sp.copy(K)
ccpsi = sp.zeros(nodes)

for i in range(psi_boundary):
    index = int(malha.dirichlet_points[i][0]-1)
    value = psi_dirichlet[index]
    for j in range(nodes):
        ccpsi[j] -= value * K[j, index]
        if j != index:
            K_psi[index, j] = 0
            K_psi[j, index] = 0
        else:
            K_psi[index, j] = 1

F_psi = sp.dot(M, Wz_old) + ccpsi
for i in range(psi_boundary):
    index = int(malha.dirichlet_points[i][0]-1)
    F_psi[index] = psi_dirichlet[index]

Psi_old = sp.linalg.solve(K_psi, F_psi)


cc_temp = sp.zeros(nodes)
len_neu = len(malha.Boundary_Neumann)
for i in range(len_neu):
    node1 = malha.Boundary_Neumann[i][0]
    node2 = malha.Boundary_Neumann[i][1]
    for j in range(len(malha.neumann_points)):
        if node1 == int(malha.neumann_points[j,0]):
            index1 = j
        if node2 == int(malha.neumann_points[j,0]):
            index2 = j
    value1 = malha.neumann_points[index1][1]
    value2 = malha.neumann_points[index2][1]
    length = sp.sqrt((x[node1]-x[node2])**2 + (y[node1]-y[node2])**2)
    neu_bc1 = (length/2)*value1
    neu_bc2 = (length/2)*value2
    cc_temp[node1-1] -= neu_bc1
    cc_temp[node2-1] -= neu_bc2

T_a = sp.zeros(nodes)
u_a = sp.zeros(nodes)
for i in range(nodes):
    T_a[i] = (-15./48. + 2.*(malha.X[i]) + 1.5*(malha.Y[i]-0.5)**2-(malha.Y[i]-0.5)**4)
    u_a[i] = 6.0*(malha.Y[i]-malha.Y[i]**2)

# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

for t in range(0, tempo-1):
    print "Solving System " + str((float(t)/(tempo-1))*100) + "%"

    for i in range(psi_boundary):
        index = int(malha.dirichlet_points[i][0]-1)
        vy[index] = 0.0
        vx[index] = 0.0

    for i in range(velx_quant):
        index = int(velx_nodes[i])
        vx[index] = 1.0

    # B.C. Vorticidade
    Wcc = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))
    ccomega = sp.zeros(nodes)

    # Solução de Wz e Psi
    #Conv = sp.dot(sp.diag(vx), Gx) + sp.dot(sp.diag(vy), Gy)
    Conv = vx* Gx + vy* Gy
    LHS_Re = Mdt + K_Re + Conv
    LHS_omega = sp.copy(LHS_Re)

    for i in range(psi_boundary):
        index = int(malha.dirichlet_points[i][0]-1)
        value = Wcc[index]
        for j in range(nodes):
            ccomega[j] -= value * LHS_Re[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    F_omega = sp.dot(Mdt, Wz_old) + ccomega #/ Re #- sp.dot(Conv, Wz_old)

    for i in range(psi_boundary):
        index = int(malha.dirichlet_points[i][0]-1)
        F_omega[index] = Wcc[index]

    Wz_new = sp.linalg.solve(LHS_omega, F_omega)


    F_psi = sp.dot(M, Wz_new) + ccpsi
    for i in range(psi_boundary):
        index = int(malha.dirichlet_points[i][0]-1)
        F_psi[index] = psi_dirichlet[index]

    Psi_new = sp.linalg.solve(K_psi, F_psi)

    # Temperatura

    LHS_T = Mdt + Conv + (1.0/(Re*Pr))*K
    LHS_temp = sp.copy(LHS_T)

    cctemp = sp.zeros(nodes)
    for i in range(velx_quant):
        index = int(velx_nodes[i])
        value = 0.0
        for j in range(nodes):
            cctemp[j] -= value * LHS_T[j, index]
            if j != index:
                LHS_temp[index, j] = 0
                LHS_temp[j, index] = 0
            else:
                LHS_temp[index, j] = 1

    F_temp = sp.dot(Mdt, Temp_old) + cctemp + cc_temp / (Re*Pr)

    for i in range(velx_quant):
        index = int(velx_nodes[i])
        F_temp[index] = 0.0

    Temp_new = sp.linalg.solve(LHS_temp, F_temp)

    # Salvar VTK
    vtk = Io.InOut(x, y, ien, len(x), len(ien), Temp_old, Psi_old, Wz_old, T_a, u_a, vx, vy)
    vtk.saveVTK(cwd+"/results", arquivo + str(t+1))

    Psi_old = sp.copy(Psi_new)
    Wz_old = sp.copy(Wz_new)
    Temp_old = sp.copy(Temp_new)

    # Calculo de Vx e Vy
    vx = sp.dot(MinvLump, sp.dot(Gy, Psi_new))
    vy = -1.0 * sp.dot(MinvLump, sp.dot(Gx, Psi_new))

#----------------- Fim de Loop -------------------
#-------------------------------------------------
