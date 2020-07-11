# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
import os
import Elements as ele
import semiLagrangean as sl

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

mesh_file = "conjugate"

savename = mesh_file

fluid_mesh = gm.GMesh("mesh/" + mesh_file + "-fld.msh")
global_mesh = gm.GMesh("mesh/" + mesh_file + ".msh")

x_fluid = fluid_mesh.X
y_fluid = fluid_mesh.Y + 1
ien_fluid = fluid_mesh.IEN
nodes_fluid = len(x_fluid)
num_ele_fluid = len(ien_fluid)

axisym_tri = ele.Linear(x_fluid, y_fluid)

x_global = global_mesh.X
y_global = global_mesh.Y + 1
ien_global = global_mesh.IEN
nodes_global = len(x_global)
num_ele_global = len(ien_global)

axisym_tri_global = ele.Linear(x_global, y_global)

count = 0
Vel_points_convert = sp.zeros(nodes_global) - 1
for i in range(nodes_global):
    for j in range(nodes_fluid):
        if (x_global[i] == x_fluid[j]) and (y_global[i] == y_fluid[j]):
            Vel_points_convert[i] = j
            count += 1
            print count, " / ", nodes_fluid


# -------------------------------------------------------
#     Simulation Parameters
# -------------------------------------------------------

dt = 0.05
time = 1000

# Fluid properties
rho_fluid = 1000
viscosity_din = 1
viscosity_kin = 0.001
alpha_fluid = 0.001
alpha_solid = 0.05

# -------------------------------------------------------
#     Initial Variables
# -------------------------------------------------------

# Flow
vz = sp.zeros(nodes_fluid)
vr = sp.zeros(nodes_fluid)

# -------------------------------------------------------
#     Matrix Assembly
# -------------------------------------------------------


def fem_matrix(_x, _y, _numele, _numnode, _ien, _alfa, _axi):
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr2_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m1r_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gxr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx = sp.zeros((_numnode, _numnode), dtype="float64")
    gy = sp.zeros((_numnode, _numnode), dtype="float64")
    gyr_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):

        alfa = 0
        for i in range(3):
            alfa += _alfa[_ien[elem, i]]

        alfa = alfa * (1/3.0)

        v = [_ien[elem, 0], _ien[elem, 1], _ien[elem, 2]]
        _axi.getAxiSym(v)
        ele_radius = (_y[v[0]] + _y[v[1]] + _y[v[2]])/3.

        for i_local in range(3):
            i_global = _ien[elem, i_local]

            for j_local in range(3):
                j_global = _ien[elem, j_local]

                k_global[i_global, j_global] += alfa * ele_radius * axisym_tri.kxx[i_local, j_local] + \
                                                alfa * ele_radius * axisym_tri.kyy[i_local, j_local]

                mr_global[i_global, j_global] += ele_radius*axisym_tri.mass[i_local, j_local]
                mr2_global[i_global, j_global] += (ele_radius**2)*axisym_tri.mass[i_local, j_local]
                m1r_global[i_global, j_global] += (1./ele_radius)*axisym_tri.mass[i_local, j_local]
                m_global[i_global, j_global] += axisym_tri.mass[i_local, j_local]
                gxr_global[i_global, j_global] += ele_radius*axisym_tri.gx[i_local, j_local]
                gx[i_global, j_global] += axisym_tri.gx[i_local, j_local]
                gyr_global[i_global, j_global] += ele_radius*axisym_tri.gy[i_local, j_local]
                gy[i_global, j_global] += axisym_tri.gy[i_local, j_local]

    return k_global, m_global, mr_global, mr2_global, m1r_global, gxr_global, gyr_global, gx, gy

alfa_global = sp.zeros(nodes_global)
for i in range(nodes_global):
    if y_global[i] > min(y_fluid):
        alfa_global[i] = alpha_fluid
    else:
        alfa_global[i] = alpha_solid

print " Fluid Matrix Assembly"
K, M, Mr, Mr2, M1r, Gxr, Gyr, Gx, Gy = fem_matrix(x_fluid, y_fluid, num_ele_fluid,
                                                  nodes_fluid, ien_fluid, sp.ones(nodes_fluid), axisym_tri)

print " Global Matrix Assembly"

Kh, Mh, Mrh, Mr2h, M1rh, Gxrh, Gyrh, Gxh, Gyh = fem_matrix(x_global, y_global, num_ele_global,
                                                           nodes_global, ien_global, alfa_global, axisym_tri_global)

Mdth = Mrh/dt
Mdt = Mr/dt

print " Lumped Fluid"
MLumpr = sp.zeros(nodes_fluid)
MinvLumpr = sp.zeros(nodes_fluid)
for i in range(nodes_fluid):
    for j in range(nodes_fluid):
        MLumpr[i] += Mr[i, j]
    MinvLumpr[i] = 1. / MLumpr[i]

print " Lumped Global"

temp_last = sp.zeros(nodes_global)

MLumprh = sp.zeros(nodes_global)
MinvLumprh = sp.zeros(nodes_global)
for i in range(nodes_global):
    if y_global[i]<= min(y_fluid):
        temp_last[i] = 10
    else:
        temp_last[i] = 1
    for j in range(nodes_global):
        MLumprh[i] += Mrh[i, j]
    MinvLumprh[i] = 1. / MLumprh[i]

# ---------------------------------------
# Boundary and Initial Condition
# ---------------------------------------

dirichlet_len_fluid = len(fluid_mesh.dirichlet_points)

Fluid_Boundary_in = sp.zeros(nodes_fluid)
Fluid_Boundary_out = sp.zeros(nodes_fluid)
Fluid_Boundary_wall = sp.zeros(nodes_fluid)
Fluid_Boundary_axis = sp.zeros(nodes_fluid)

Fluid_Boundary = sp.zeros(nodes_fluid)

list1 = []
list2 = []
list3 = []
list4 = []
list5 = []

for i in range(nodes_fluid):
    if x_fluid[i] == 0:
        Fluid_Boundary_in[i] = i
    else:
        list1.append(i)

    if x_fluid[i] == 5:
        Fluid_Boundary_out[i] = i
    else:
        list2.append(i)

    if y_fluid[i] == max(y_fluid):
        Fluid_Boundary_wall[i] = i
    else:
        list3.append(i)

    if y_fluid[i] == min(y_fluid):
        Fluid_Boundary_axis[i] = i
    else:
        list4.append(i)

    if x_fluid[i] == 0 or x_fluid[i] == 5 \
            or y_fluid[i] == min(y_fluid) or y_fluid[i] == max(y_fluid):
        Fluid_Boundary[i] = i
    else:
        list5.append(i)


Fluid_Boundary_in = sp.delete(Fluid_Boundary_in, list1, axis=0)
Fluid_Boundary_out = sp.delete(Fluid_Boundary_out, list2, axis=0)
Fluid_Boundary_wall = sp.delete(Fluid_Boundary_wall, list3, axis=0)
Fluid_Boundary_axis = sp.delete(Fluid_Boundary_axis, list4, axis=0)

Fluid_Boundary = sp.delete(Fluid_Boundary, list5, axis=0)


num_omega = len(Fluid_Boundary)


vz_a = sp.zeros(nodes_fluid)
psi_a = sp.zeros(nodes_fluid)
omega_a = sp.zeros(nodes_fluid)
dpdx = -16
gg = -dpdx/(4 * viscosity_din)

# Inside cylinder
# for i in range(nodes_fluid):
#     vz_a[i] = gg * (1 - y_fluid[i]**2)
#     psi_a[i] = gg * (0.5 * y_fluid[i]**2 - 0.25 * y_fluid[i]**4)
#     omega_a[i] = 2 * gg * y_fluid[i]

# Between cylinders
R1 = min(y_fluid)
R2 = max(y_fluid)
for i in range(nodes_fluid):
    ri = y_fluid[i]
    vz_a[i] = gg * (R1**2 - ri ** 2) + \
              gg * (R2**2 - R1**2) * (sp.log(ri/R1) / sp.log(R2/R1))

    c0 = gg * 0.25 * R1**4 - gg * (R2**2 - R1**2) * (1./sp.log(R2/R1)) * \
        0.5 * R1**2 * 0.5
    psi_a[i] = gg * (0.5 * R1**2 * ri ** 2 - 0.25 * ri**4) + \
        gg * (R2**2 - R1**2) * (1./sp.log(R2/R1)) * 0.5 * ri**2 * \
        (sp.log(ri) - 0.5 - sp.log(R1)) - c0

    omega_a[i] = 2 * gg * ri - gg * (R2**2 - R1**2) * (1./sp.log(R2/R1)) * (1./ri)


# --------------------------------------
# Psi K matrix with Dirichlet BC -- K_psi

print " LHS_psi"


K_psi = K + 2 * Gy
ccpsi = sp.zeros(nodes_fluid)
for i in range(dirichlet_len_fluid):
    index = int(fluid_mesh.dirichlet_points[i][0] - 1)
    # value = fluid_mesh.dirichlet_points[i][1]               # Inside cylinder
    value = psi_a[index]       # Between Cylinders
    for j in range(nodes_fluid):
        ccpsi[j] -= value * K[j, index]
        if j != index:
            K_psi[index, j] = 0
            K_psi[j, index] = 0
        else:
            K_psi[index, j] = 1
# --------------------------------------


for i in Fluid_Boundary:
    j = int(i)
    vr[j] = 0.0
    if x_fluid[j] == 0 and min(y_fluid) < y_fluid[j] < max(y_fluid):
        vz[j] = 2.0
        # vz[j] = vz_a[j]
    if y_fluid[j] == 1:
        vz[j] = 0
        vr[j] = 0


omega_last = sp.multiply(MinvLumpr, (sp.dot(Gxr, vr) - sp.dot(Gyr, vz)))

F_psi = sp.dot(Mr2, omega_last) + ccpsi
for i in range(dirichlet_len_fluid):
    index = int(fluid_mesh.dirichlet_points[i][0] - 1)
    # F_psi[index] = fluid_mesh.dirichlet_points[i][1]
    F_psi[index] = psi_a[index]

psi = linalg.solve(K_psi, F_psi)

# --------------------------------------
# K temperature

print " LHS_t"

dirichlet_len_global = len(global_mesh.dirichlet_points)


LHS_temp = Mdth + Kh
Var = sp.copy(LHS_temp)
cctemp = sp.zeros(nodes_global)
for i in range(dirichlet_len_global):
    index = int(global_mesh.dirichlet_points[i][0] - 1)
    value = global_mesh.dirichlet_points[i][1]
    for j in range(nodes_global):
        cctemp[j] -= value * Var[j, index]
        if j != index:
            LHS_temp[index, j] = 0
            LHS_temp[j, index] = 0
        else:
            LHS_temp[index, j] = 1

vz_global = sp.zeros(nodes_global)
vr_global = sp.zeros(nodes_global)
# --------------------------------------


# ----------------------------------------------------------
# -------------------- Time iteration ----------------------

print " neighbours"

neighbour = sl.neighbourElements(nodes_fluid, ien_fluid)
neighbour_global = sl.neighbourElements(nodes_global, ien_global)
for t in range(0, time):
    print "Solving System " + str((float(t)/time)*100) + "%"

    for i in range(num_omega):
        index = int(Fluid_Boundary[i])
        vr[index] = 0
        if y_fluid[index] == max(y_fluid):
            vz[index] = 0.
            vr[index] = 0.
        if x_fluid[index] == 0 and min(y_fluid) < y_fluid[index] < max(y_fluid):
            vz[index] = 2.
            # vz[index] = vz_a[index]


    for i in range(nodes_global):
        if Vel_points_convert[i] >= 0:
            vz_global[i] = vz[int(Vel_points_convert[i])]
            vr_global[i] = vr[int(Vel_points_convert[i])]
        else:
            vz_global[i] = 0
            vr_global[i] = 0

    omega_dep = sl.Linear2D(nodes_fluid, neighbour, ien_fluid, x_fluid,
                            y_fluid, vz, vr, dt, omega_last)

    # B.C. Vorticidade
    W_in = sp.multiply(MinvLumpr, (sp.dot(Gxr, vr) - sp.dot(Gyr, vz)))
    W_axis = W_in
    W_wall = W_in
    ccomega = sp.zeros(nodes_fluid)

    # Solução de Wz e Psi
    Conv = vr * sp.diag(Mr)

    LHS_Ni = Mdt + K * viscosity_kin - sp.diag(Conv) + M1r * viscosity_kin
    LHS_omega = sp.copy(LHS_Ni)

    for i in range(len(Fluid_Boundary_in)):
        index = int(Fluid_Boundary_in[i])
        # value = 8*y_fluid[index]
        value = omega_a[index]
        for j in range(nodes_fluid):
            ccomega[j] -= value * LHS_Ni[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    for i in range(len(Fluid_Boundary_axis)):
        index = int(Fluid_Boundary_axis[i])
        # value = W_axis[index]
        value = omega_a[index]
        for j in range(nodes_fluid):
            ccomega[j] -= value * LHS_Ni[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    for i in range(len(Fluid_Boundary_wall)):
        index = int(Fluid_Boundary_wall[i])
        # value = W_wall[index]
        value = omega_a[index]
        for j in range(nodes_fluid):
            ccomega[j] -= value * LHS_Ni[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    F_omega = sp.dot(Mdt, omega_dep) + ccomega

    for i in range(len(Fluid_Boundary_in)):
        index = int(Fluid_Boundary_in[i])
        # F_omega[index] = W_in[index]
        F_omega[index] = omega_a[index]


    for i in range(len(Fluid_Boundary_wall)):
        index = int(Fluid_Boundary_wall[i])
        # F_omega[index] = W_wall[index]
        F_omega[index] = omega_a[index]

    for i in range(len(Fluid_Boundary_axis)):
        index = int(Fluid_Boundary_axis[i])
        # F_omega[index] = W_axis[index]
        F_omega[index] = omega_a[index]

    omega = linalg.solve(LHS_omega, F_omega)

    F_psi = sp.dot(Mr2, omega) + ccpsi
    for i in range(dirichlet_len_fluid):
        index = int(fluid_mesh.dirichlet_points[i][0]-1)
        # F_psi[index] = fluid_mesh.dirichlet_points[i][1]
        F_psi[index] = psi_a[index]

    psi = linalg.solve(K_psi, F_psi)

    temp_dep = sl.Linear2D(nodes_global, neighbour_global, ien_global, x_global,
                           y_global, vz_global, vr_global, dt, temp_last)

    # --------------------------------------
    # RHS matrix with Dirichlet BC

    RHS_temp = sp.dot(Mdth, temp_dep) + cctemp
    for i in range(dirichlet_len_global):
        index = int(global_mesh.dirichlet_points[i][0] - 1)
        RHS_temp[index] = global_mesh.dirichlet_points[i][1]

    temp = linalg.solve(LHS_temp, RHS_temp)


    # Salvar VTK
    vtk_t = IO.InOut(x_global, y_global, ien_global, nodes_global, num_ele_global, temp, None, None
                     , None, None, vz_global, vr_global)
    vtk_t.saveVTK(cwd+"/results", savename + str(t+1))

    vtk = IO.InOut(x_fluid, y_fluid, ien_fluid, nodes_fluid, num_ele_fluid, psi, omega, vz_a
                     , psi_a, omega_a, vz, vr)
    vtk.saveVTK(cwd+"/results", savename + "-fld" + str(t+1))

    omega_last = sp.copy(omega)
    temp_last = sp.copy(temp

                        )
    # Calculo de vz e vr

    vz = sp.multiply(MinvLumpr, sp.dot(Gy, psi))
    vr = -1.0 * sp.multiply(MinvLumpr, sp.dot(Gx, psi))


# ----------------- Fim de Loop -------------------
# -------------------------------------------------