# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import os
import semiLagrangean as sl

cwd = os.getcwd()

arquivo = "plane"

malha = Gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN
nodes = len(x)
num_ele = len(ien)

dt = 0.01
tempo = 100
Re = 10

<<<<<<< HEAD
p_lagrange = 0.2
p_smooth = 1.
=======
p_lagrange = 0.4
p_smooth = 0.8
>>>>>>> e1d3b68bdae3f65dfaf0bfa5b349e426c22a65a4
p_wave = 0


def smoothMesh(_neighbour_nodes, _mesh, _x, _y, _dt):

    xx = sp.copy(_x)
    yy = sp.copy(_y)
    vx_disp = sp.zeros(len(xx))
    vy_disp = sp.zeros(len(xx))
    for i in range(len(_neighbour_nodes)):

        flag = 0
        for k in _mesh.Boundary_Nodes:
            if i == k:
                flag = 1
                break

        if flag == 1:
            continue

        vertex_position = sp.array([xx[i], yy[i]])
        nghN = _neighbour_nodes[i]
        num_nghb = len(nghN)
        distance_vectors = sp.zeros((num_nghb, 2))
        displacement_vector = sp.zeros(2)
        for j in range(num_nghb):
            _index = nghN[j]
            distance_vectors[j][0] = xx[_index] - vertex_position[0]
            distance_vectors[j][1] = yy[_index] - vertex_position[1]
            displacement_vector += distance_vectors[j]
        displacement_vector *= (1./num_nghb)
        displacement_velocity = displacement_vector / dt
        vx_disp[i] = displacement_velocity[0]
        vy_disp[i] = displacement_velocity[1]

      #  print nghN, "\n", distance_vectors, "\n", displacement_vector, "\n", new_position

    return vx_disp, vy_disp


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

K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

K_Re = K/Re
Mdt = M/dt
MLump = sp.zeros((nodes,nodes))

Boundary = sp.zeros(nodes)
lista = []
for i in range(nodes):
    if x[i] == 0.0 or y[i] == 0.0 or y[i] == 1.0 or x[i] == 5.0:
        Boundary[i] = i
    else:
        lista.append(i)
    for j in range(nodes):
        MLump[i, i] += M[i, j]

MinvLump = linalg.inv(MLump)

Boundary = sp.delete(Boundary, lista, axis=0)
num_bc = len(Boundary)

# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------
psi_dirichlet = sp.zeros(nodes)
bc_omega = sp.zeros(nodes)
for i in Boundary:
    j = int(i)
    vy[j] = 0.0
    if y[j] == 1.0:
        psi_dirichlet[j] = 1.0
        vx[j] = 0.0
    if y[j] == 0.0:
        psi_dirichlet[j] = 0.0
        vx[j] = 0.0
    if x[j] == 0.0 and (y[j] > 0.0 and y[j] < 1.0):
        vx[j] = 1.0

Wz_old = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))

K_psi = sp.copy(K)
ccpsi = sp.zeros(nodes)

for i in range(num_bc):
    index = int(Boundary[i])
    value = psi_dirichlet[index]
    if y[index] == 0.0 or y[index] == 1.0:
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
    if y[index] == 0.0 or y[index] == 1.0:
        F_psi[index] = psi_dirichlet[index]

Psi_old = sp.linalg.solve(K_psi, F_psi)

sp.random.seed(1)

# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------
neighbour_ele, neighbour_nodes = sl.neighbourElements2(nodes, ien)
for t in range(0, tempo-1):
    print "Solving System " + str((float(t)/(tempo-1))*100) + "%"

    vx_smooth, vy_smooth = smoothMesh(neighbour_nodes, malha, x, y, dt)
    vx_wave = sp.random.rand(nodes) * 0.2*sp.cos(2*sp.pi*5*t*dt)
    vy_wave = sp.random.rand(nodes) * 0.2*sp.cos(2*sp.pi*5*t*dt)

    vxAle = p_lagrange * vx + p_smooth * vx_smooth + p_wave * vx_wave
    vyAle = p_lagrange * vy + p_smooth * vy_smooth #+ p_wave * vy_wave

    for i in range(num_bc):
        index = int(Boundary[i])
        vy[index] = 0.0
        vxAle[index] = 0
        vyAle[index] = 0
        if y[index] == max(y) or y[index] == min(y):
            vx[index] = 0.0
        if x[index] == min(x) and max(y) > y[index] > min(y):
            vx[index] = 1.0

    vx_sl = vx - vxAle
    vy_sl = vy - vyAle
    #Wz_dep = sl.Linear2D(nodes, neighbour, ien, x, y, vx, vy, dt, Wz_old)
    Wz_dep = sl.Linear2D(nodes, neighbour_ele, ien, x, y, vx_sl, vy_sl, dt, Wz_old)


    # Solução de Wz e Psi

    x = x + vxAle * dt
    y = y + vyAle * dt
    K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)


<<<<<<< HEAD
    MinvLump = linalg.inv(MLump)
=======
    MinvLump = linalg.inv(M)
>>>>>>> e1d3b68bdae3f65dfaf0bfa5b349e426c22a65a4

    # B.C. Vorticidade
    Wcc = sp.dot(MinvLump, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))
    ccomega = sp.zeros(nodes)
    #Wcc = sp.copy(bc_omega)

    LHS = M / dt + K / Re
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

    F_omega = sp.dot(M / dt, Wz_dep) + ccomega / Re

    for i in range(num_bc):
        index = int(Boundary[i])
        F_omega[index] = Wcc[index]

    Wz_new = sp.linalg.solve(LHS_omega, F_omega)

    K_psi = sp.copy(K)
    ccpsi = sp.zeros(nodes)

    for i in range(num_bc):
        index = int(Boundary[i])
        value = psi_dirichlet[index]
        if y[index] == 0.0 or y[index] == 1.0:
            for j in range(nodes):
                ccpsi[j] -= value * K[j, index]
                if j != index:
                    K_psi[index, j] = 0
                    K_psi[j, index] = 0
                else:
                    K_psi[index, j] = 1

    F_psi = sp.dot(M, Wz_new) + ccpsi
    for i in range(num_bc):
        index = int(Boundary[i])
        if y[index] == 0.0 or y[index] == 1.0:
            F_psi[index] = psi_dirichlet[index]

    Psi_new = sp.linalg.solve(K_psi, F_psi)

    # Salvar VTK
    vtk = Io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, Wz_dep, None, None, vx, vy, vx_sl, vy_sl)
    vtk.saveVTK(cwd+"/results", arquivo + str(t+1))

    Psi_old = sp.copy(Psi_new)
    Wz_old = sp.copy(Wz_new)

    # Calculo de Vx e Vy
    vx = sp.dot(MinvLump, sp.dot(Gy, Psi_new))
    vy = -1.0 * sp.dot(MinvLump, sp.dot(Gx, Psi_new))

#----------------- Fim de Loop -------------------
#-------------------------------------------------
