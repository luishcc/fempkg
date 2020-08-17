# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
import matplotlib.pyplot as plt
import Elements as ele
import os

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

# mesh_file = "poiseuille"
# mesh_file = "fine"
mesh_file = "poi-var"

savename = mesh_file + "-test"

fluid_mesh = gm.GMesh("mesh/" + mesh_file + "-fld.msh")

x_fluid = fluid_mesh.X
y_fluid = fluid_mesh.Y
# y_fluid = fluid_mesh.Y + 1
ien_fluid = fluid_mesh.IEN
nodes_fluid = len(x_fluid)
num_ele_fluid = len(ien_fluid)

axisym_tri = ele.Linear(x_fluid, y_fluid)



# -------------------------------------------------------
#     Simulation Parameters
# -------------------------------------------------------

# Fluid properties
rho_fluid = 1000
viscosity_din = 1
viscosity_kin = 1

# -------------------------------------------------------
#     Initial Variables
# -------------------------------------------------------

# Flow
psi_numeric = sp.zeros(nodes_fluid)
vz = sp.zeros(nodes_fluid)
vr = sp.zeros(nodes_fluid)

# -------------------------------------------------------
#     Matrix Assembly
# -------------------------------------------------------


def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    k2_global = sp.zeros((_numnode, _numnode), dtype="float64")
    kpsi_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr2_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gxr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx = sp.zeros((_numnode, _numnode), dtype="float64")
    gy = sp.zeros((_numnode, _numnode), dtype="float64")
    gyr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    xx = sp.zeros(_numnode, dtype="float64")
    yy = sp.zeros(_numnode, dtype="float64")


    for elem in range(_numele):

        v = [_ien[elem, 0], _ien[elem, 1], _ien[elem, 2]]
        axisym_tri.getAxiSym(v)

        ele_radius = (_y[v[0]] + _y[v[1]] + _y[v[2]])/3.

        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]

        a1 = xx[0] * yy[1] - xx[1] * yy[0]
        a2 = xx[2] * yy[0] - xx[0] * yy[2]
        a3 = xx[1] * yy[2] - xx[2] * yy[1]
        area = (a1 + a2 + a3) / 2.

        if area == 0:
            print elem

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]

                k_global[i_global, j_global]     += ele_radius*axisym_tri.kxx[i_local, j_local]+\
                                                    ele_radius*axisym_tri.kyy[i_local, j_local]+\
                                                    axisym_tri.gy[i_local, j_local]

                k2_global[i_global, j_global] += ele_radius * axisym_tri.kxx[i_local, j_local] + \
                                                 ele_radius * axisym_tri.kyy[i_local, j_local] + \
                                                 2.*axisym_tri.gy[i_local, j_local]

                #k2_global[i_global, j_global] += axisym_tri.kxx[i_local, j_local] + \
                 #                                axisym_tri.kyy[i_local, j_local]+ \
                  #                              2*axisym_tri.gy[i_local, j_local]/ele_radius


                mr_global[i_global, j_global]    += ele_radius*axisym_tri.mass[i_local, j_local]
                #mr_global[i_global, j_global]    += axisym_tri.mr[i_local, j_local]
                #mr_global[i_global, j_global]    += ele_radius*m_local[i_local, j_local]*(area/6.)
                mr2_global[i_global, j_global]   += (ele_radius**2)*axisym_tri.mass[i_local, j_local]
                m_global[i_global, j_global]     += axisym_tri.mass[i_local, j_local]
                #m_global[i_global, j_global]     += m_local[i_local, j_local] * (area/6.)
                gxr_global[i_global, j_global]   += ele_radius*axisym_tri.gx[i_local, j_local]
                gx[i_global, j_global]           += axisym_tri.gx[i_local, j_local]
                gyr_global[i_global, j_global]   += ele_radius*axisym_tri.gy[i_local, j_local]
                gy[i_global, j_global]           += axisym_tri.gy[i_local, j_local]

    return k2_global, m_global, mr_global, gxr_global, gyr_global, gx, gy, mr2_global

print "Matrix Assembly"
K, M, Mr, Gxr, Gyr, Gx, Gy, Mr2 = fem_matrix(x_fluid, y_fluid, num_ele_fluid, nodes_fluid, ien_fluid)

print "Lumped Mass"
MLump = sp.zeros(nodes_fluid)
MinvLump = sp.zeros(nodes_fluid)
for i in range(nodes_fluid):
    for j in range(nodes_fluid):
        MLump[i] += Mr[i, j]
    MinvLump[i] = 1. / MLump[i]

# ---------------------------------------
# Boundary and Initial Condition
# ---------------------------------------

dirichlet_len_fluid = len(fluid_mesh.dirichlet_points)

vz_a = sp.zeros(nodes_fluid)
psi_a = sp.zeros(nodes_fluid)
psi2 = sp.zeros(nodes_fluid)
omega_a = sp.zeros(nodes_fluid)
dpdx = -16
gg = -dpdx/(4 * viscosity_din)

# Inside cylinder
for i in range(nodes_fluid):
    vz_a[i] = gg * (1 - y_fluid[i]**2)
    psi_a[i] = gg * (0.5 * y_fluid[i]**2 - 0.25 * y_fluid[i]**4)
    omega_a[i] = 2 * gg * y_fluid[i]

# Between cylinders
# R1 = min(y_fluid)
# R2 = max(y_fluid)
# for i in range(nodes_fluid):
#     ri = y_fluid[i]
#     vz_a[i] = gg * (R1 - ri ** 2) + \
#               gg * (R2**2 - R1**2) * (sp.log(ri/R1) / sp.log(R2/R1))
#
#     c0 = gg * 0.25 * R1**4 - gg * (R2**2 - R1**2) * (1./sp.log(R2/R1)) * \
#          0.5 * R1**2 * 0.5
#     psi_a[i] = gg * (0.5 * R1**2 * ri ** 2 - 0.25 * ri**4) + \
#                gg * (R2**2 - R1**2) * (1./sp.log(R2/R1)) * 0.5 * ri**2 * \
#                (sp.log(ri) - 0.5 - sp.log(R1)) - c0
#
#     omega_a[i] = 2 * gg * ri - gg * (R2**2 - R1**2) * (1./sp.log(R2/R1)) * (1./ri)

# --------------------------------------
# Psi K matrix with Dirichlet BC -- K_psi

print "LHS Psi eq."

#K_psi = K + 2 * Gy
K_psi = 1 * K
ccpsi = sp.zeros(nodes_fluid)
for i in range(dirichlet_len_fluid):
    index = int(fluid_mesh.dirichlet_points[i][0] - 1)
    value = fluid_mesh.dirichlet_points[i][1]
    # value = fluid_mesh.dirichlet_points[i][1] * 2
    for j in range(nodes_fluid):
        ccpsi[j] -= value * K[j, index]
        if j != index:
            K_psi[index, j] = 0
            K_psi[j, index] = 0
        else:
            K_psi[index, j] = 1

F_psi = sp.dot(Mr2, omega_a) + ccpsi
for i in range(dirichlet_len_fluid):
    index = int(fluid_mesh.dirichlet_points[i][0] - 1)
    F_psi[index] = fluid_mesh.dirichlet_points[i][1]
    # F_psi[index] = fluid_mesh.dirichlet_points[i][1]*2


psi_numeric = sp.linalg.solve(K_psi, F_psi)

# Calculo de vz e vr

vz = sp.multiply(MinvLump, sp.dot(Gy, psi_numeric))
# vz = sp.multiply(MinvLump, sp.dot(Gy, psi_numeric))

vz_axis = sp.multiply(MinvLump, sp.dot(K, psi_numeric))
# vz_axis = vz


vr = -1 * sp.multiply(MinvLump, sp.dot(Gx, psi_numeric))

for i in range(nodes_fluid):
    if y_fluid[i] == 0:
        print i
        vz[i] = vz_axis[i]

    if y_fluid[i] == max(y_fluid):
        print i
        vz[i] = 0


#omega = -sp.multiply(MinvLump,  sp.dot(Gyr, vz)) + sp.multiply(MinvLump,  sp.dot(Gxr, vr))
omega = -sp.multiply(MinvLump,  sp.dot(Gyr, vz_a))
#omega = sp.linalg.solve(Mr, -sp.dot(Gyr, vz_a))

erro = abs(vz_a - vz)
erro_psi = abs(psi_a - psi_numeric)

# Salvar VTK
vtk_t = IO.InOut(x_fluid, y_fluid, ien_fluid, nodes_fluid, num_ele_fluid, psi_numeric, omega, vz_a
                 , psi_a, omega_a, vz, vr)
vtk_t.saveVTK(cwd+"/results", savename + "Test")

