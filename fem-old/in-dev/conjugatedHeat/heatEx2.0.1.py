# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import os

cwd = os.getcwd()


# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

arquivo = "micro"

malha_total = Gm.GMesh(arquivo+".msh")

x_total = malha_total.X
y_total = malha_total.Y
ien_total = malha_total.IEN
nodes_total = len(x_total)
num_ele_total = len(ien_total)

malha_fldH = Gm.GMesh(arquivo+"-fldh.msh")

x_fldH = malha_fldH.X
y_fldH = malha_fldH.Y
ien_fldH = malha_fldH.IEN
nodes_fldH = len(x_fldH)
num_ele_fldH = len(ien_fldH)

for i in range(len(malha_total.dirichlet_points)):
    point = int(malha_total.dirichlet_points[i][0]-1)
    if (x_total[point] == 0) and (0.18 < y_total[point] < 0.21):
        malha_total.dirichlet_points[i][1] = 60

Vel_points_convertH = sp.zeros(nodes_total) -1
for i in range(nodes_total):
    for j in range(nodes_fldH):
        if (x_total[i] == x_fldH[j]) and (y_total[i] == y_fldH[j]):
            Vel_points_convertH[i] = j


# -------------------------------------------------------
#     Defining Parameters
# -------------------------------------------------------

dt = 0.5
tempo = 500

Ratio_k = 5

Ni_H = 0.00015
TermCond_fld = 10
Rho_fld = 1000.0
Cp_fld = 25.0
alfa_fldH = TermCond_fld / (Rho_fld * Cp_fld)

Rho_solid = 1000.0
Cp_solid = 25.0
TermCond_solid = Ratio_k * TermCond_fld

alfa_solid = TermCond_solid / (Cp_solid * Rho_solid)


alfa_solid = 1.0
alfa_fldH = 0.001
alfa_solid2 = 1.2



# ---------------------------------------
# Wz, Psi, T e velocidade inicial
# ---------------------------------------

Psi_new_h = sp.zeros(nodes_fldH, dtype="float64")
Wz_new_h = sp.zeros(nodes_fldH, dtype="float64")
vx_h = sp.zeros(nodes_fldH, dtype="float64")
vy_h = sp.zeros(nodes_fldH, dtype="float64")

vx_total = sp.zeros(nodes_total, dtype="float64")
vy_total = sp.zeros(nodes_total, dtype="float64")

Temp_old = sp.ones(nodes_total, dtype="float64")*20
for i in range(len(malha_total.dirichlet_points)):
    index = int(malha_total.dirichlet_points[i][0]-1)
    value = malha_total.dirichlet_points[i][1]
    Temp_old[index] = value


# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien, _alfa):
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
                gx_local[i,j] = b[j] * (1/6.)
                gy_local[i,j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]*alfa
                m_global[i_global, j_global] += m_local[i_local, j_local]* (Area/12.)
                gx_global[i_global, j_global] += gx_local[i_local, j_local]
                gy_global[i_global, j_global] += gy_local[i_local, j_local]


    return  k_global, m_global, gx_global, gy_global

# --------------------------------------
# Total Matrices (Heat equation)


Alfa_total = sp.zeros(nodes_total)
for i in range(nodes_total):
    if y_total[i] > 0.2 or y_total[i] < 0.4:
        Alfa_total[i] = alfa_fldH
    if y_total[i] <= 0.2:
        Alfa_total[i] = alfa_solid
    if y_total[i] >= 0.4:
        Alfa_total[i] = alfa_solid2
    if ((0.125 < x_total[i] < 0.25) or (0.375 < x_total[i] < 0.5) or (0.625 < x_total[i] < 0.75) or
            (0.875 < x_total[i] < 1.0)) and (y_total[i] > 0.38):
        Alfa_total[i] = alfa_solid2


K, M, Gx, Gy = fem_matrix(x_total, y_total, num_ele_total, nodes_total, ien_total, Alfa_total)

Mdt = M/dt
# --------------------------------------


# --------------------------------------
# Fluid Region (Psi, Omega)

Alfa_fluid_h = sp.ones(nodes_fldH)

K_fld_h, M_fld_h, Gx_fld_h, Gy_fld_h = fem_matrix(x_fldH, y_fldH, num_ele_fldH, nodes_fldH, ien_fldH, Alfa_fluid_h)

K_fld_ni_h = K_fld_h * Ni_H
Mdt_fld_h = M_fld_h/dt

MLump_h = sp.zeros((nodes_fldH, nodes_fldH))
for i in range(nodes_fldH):
    for j in range(nodes_fldH):
        MLump_h[i, i] += M_fld_h[i, j]
MinvLump_h = linalg.inv(MLump_h)

K_psi_h = sp.copy(K_fld_h)
ccpsi_h = sp.zeros(nodes_fldH)
# --------------------------------------

# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------

dirichlet_len_total = len(malha_total.dirichlet_points)
dirichlet_len_fluid_h = len(malha_fldH.dirichlet_points)
neumann_len_fluid_h = len(malha_fldH.neumann_points)

Fluid_Boundary_h = malha_fldH.dirichlet_points[:, 0]

lista = sp.zeros(neumann_len_fluid_h)
for i in range(neumann_len_fluid_h):
    lista[i] = malha_fldH.neumann_points[i][0]

Fluid_Boundary_h = sp.append(Fluid_Boundary_h, lista, axis=0) - 1

num_omega_h = len(Fluid_Boundary_h)

for i in range(num_omega_h):
    Fluid_Boundary_h[i] = int(Fluid_Boundary_h[i])


# --------------------------------------
# Psi K matrix with Dirichlet BC -- K_psi

for i in range(dirichlet_len_fluid_h):
    index = int(malha_fldH.dirichlet_points[i][0] - 1)
    value = malha_fldH.dirichlet_points[i][1]
    for j in range(nodes_fldH):
        ccpsi_h[j] -= value * K_fld_h[j, index]
        if j != index:
            K_psi_h[index, j] = 0
            K_psi_h[j, index] = 0
        else:
            K_psi_h[index, j] = 1
# --------------------------------------


# --------------------------------------
# Velocity

for i in range(nodes_fldH):
    vy_h[i] = 0.0
    vx_h[i] = 0.0

for i in range(neumann_len_fluid_h):
    index = int(malha_fldH.neumann_points[i][0]-1)
    if int(malha_fldH.neumann_points[i][2]) == 3:
        vx_h[index] = 0.1
    # if int(malha_fldH.neumann_points[i][2]) == 2:
    #     vx_h[index] = -0.1

for i in range(dirichlet_len_fluid_h):
    index = int(malha_fldH.dirichlet_points[i][0]-1)
    vx_h[index] = 0.0


# --------------------------------------

Wz_old_h = sp.dot(MinvLump_h, (sp.dot(Gx_fld_h, vy_h) - sp.dot(Gy_fld_h, vx_h)))

F_psi_h = sp.dot(M_fld_h, Wz_old_h) + ccpsi_h
for i in range(dirichlet_len_fluid_h):
    index = int(malha_fldH.dirichlet_points[i][0] - 1)
    F_psi_h[index] = malha_fldH.dirichlet_points[i][1]

Psi_old_h = sp.linalg.solve(K_psi_h, F_psi_h)


# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

check = 0
checkTemp = 0
vel_error = sp.zeros(nodes_fldH)
for t in range(0, tempo):
    print
    print "Solving System " + str((float(t)/(tempo))*100) + "%"
    print "Vorticity Error: ", check, "/", nodes_fldH*0.99
    print "Temperature Error: ", checkTemp, "/", nodes_total*0.99


    if checkTemp >= nodes_total*0.99:
        print "Converged to error in: ", t, "time steps with dt = ", dt
        exit()

    vx_h_old = sp.copy(vx_h)

    for i in range(neumann_len_fluid_h):
        index = int(malha_fldH.neumann_points[i][0]-1)
        if int(malha_fldH.neumann_points[i][2]) == 3:
            vx_h[index] = 0.1
        # if int(malha_fldH.neumann_points[i][2]) == 2:
        #     vx_h[index] = -0.1

    for i in range(dirichlet_len_fluid_h):
        index = int(malha_fldH.dirichlet_points[i][0] - 1)
        vx_h[index] = 0.0
        vy_h[index] = 0.0

    for i in range(nodes_total):
        vx_total[i] = 0
        vy_total[i] = 0
        if Vel_points_convertH[i] >= 0:
            vx_total[i] = vx_h[int(Vel_points_convertH[i])]
            vy_total[i] = vy_h[int(Vel_points_convertH[i])]

    # B.C. Vorticidade
    Wcc_h = sp.dot(MinvLump_h, (sp.dot(Gx_fld_h, vy_h) - sp.dot(Gy_fld_h, vx_h)))
    ccomega_h = sp.zeros(nodes_fldH)

    # Solução de Wz e Psi
    Conv_h = vx_h * Gx_fld_h + vy_h * Gy_fld_h
    LHS_Ni_h = Mdt_fld_h + K_fld_ni_h + Conv_h
    LHS_omega_h = sp.copy(LHS_Ni_h)

    for i in range(num_omega_h):
        index = int(Fluid_Boundary_h[i])
        value = Wcc_h[index]
        for j in range(nodes_fldH):
            ccomega_h[j] -= value * LHS_Ni_h[j, index]
            if j != index:
                LHS_omega_h[index, j] = 0
                LHS_omega_h[j, index] = 0
            else:
                LHS_omega_h[index, j] = 1

    F_omega_h = sp.dot(Mdt_fld_h, Wz_old_h) + ccomega_h  # - sp.dot(Conv_h, Wz_old_h)

    for i in range(num_omega_h):
        index = int(Fluid_Boundary_h[i])
        F_omega_h[index] = Wcc_h[index]

    if check <= nodes_fldH*0.99:
        print "solving vorticity"
        Wz_new_h = sp.linalg.solve(LHS_omega_h, F_omega_h)

    F_psi_h = sp.dot(M_fld_h, Wz_new_h) + ccpsi_h
    for i in range(dirichlet_len_fluid_h):
        index = int(malha_fldH.dirichlet_points[i][0]-1)
        F_psi_h[index] = malha_fldH.dirichlet_points[i][1]

    if check <= nodes_fldH*0.99:
        print "solving Stream Function"
        Psi_new_h = sp.linalg.solve(K_psi_h, F_psi_h)

    # Temperature
    print "calculating convection term"
    Conv_total = sp.dot(sp.diag(vx_total), Gx) + sp.dot(sp.diag(vy_total), Gy)
    #Conv_total = vx_total * Gx + vy_total * Gy

    print "calculating Heat LHS term and bc"
    LHS_T = Mdt + K + Conv_total
    LHS_temp = sp.copy(LHS_T)
    cctemp = sp.zeros(nodes_total)
    for i in range(dirichlet_len_total):
        index = int(malha_total.dirichlet_points[i][0]) - 1
        value = malha_total.dirichlet_points[i][1]
        for j in range(nodes_total):
            cctemp[j] -= value * LHS_T[j, index]
            if j != index:
                LHS_temp[index, j] = 0
                LHS_temp[j, index] = 0
            else:
                LHS_temp[index, j] = 1

    F_temp = sp.dot(Mdt, Temp_old) + cctemp
    for i in range(dirichlet_len_total):
        index = int(malha_total.dirichlet_points[i][0]) - 1
        F_temp[index] = malha_total.dirichlet_points[i][1]

    print "solving Heat"
    Temp_new = sp.linalg.solve(LHS_temp, F_temp)

    # Salvar VTK
    vtk_t = Io.InOut(x_total, y_total, ien_total, nodes_total, num_ele_total, Temp_old, Alfa_total, None, None, None, vx_total, vy_total)
    vtk_t.saveVTK(cwd+"/results", arquivo + str(t+1))

    Psi_old_h = sp.copy(Psi_new_h)
    Wz_old_h = sp.copy(Wz_new_h)


    # Calculo de Vx e Vy
    vx_h = sp.dot(MinvLump_h, sp.dot(Gy_fld_h, Psi_new_h))
    vy_h = -1.0 * sp.dot(MinvLump_h, sp.dot(Gx_fld_h, Psi_new_h))

    check = 0
    for i in range(nodes_fldH):
        if abs(vx_h[i] - vx_h_old[i]) <= 0.0000001:
            check += 1

    checkTemp = 0
    for i in range(nodes_total):
        if abs(Temp_old[i]-Temp_new[i]) <= 0.0000001:
            checkTemp += 1


    Temp_old = sp.copy(Temp_new)



# ----------------- Fim de Loop -------------------
# -------------------------------------------------
