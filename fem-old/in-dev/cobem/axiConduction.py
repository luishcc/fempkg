# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
from scipy import sparse
from scipy.sparse.linalg import spsolve
import os
import sys
import Elements as ele

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

# mesh_file = "tube"
mesh_file = "tube-coarse"
# mesh_file = "tube-coarse-dz"

# savename = mesh_file
# savename = "axi_cond_2d"
savename = "generation"

global_mesh = gm.GMesh("mesh/" + mesh_file + ".msh")

x_global = global_mesh.X
y_global = global_mesh.Y
ien_global = global_mesh.IEN
nodes_global = len(x_global)
num_ele_global = len(ien_global)

axisym_tri = ele.Linear(x_global, y_global)


# -------------------------------------------------------
#     Simulation Parameters
# -------------------------------------------------------

dt = 0.01
dt_inv = 1. / dt
time = 1000

# Fluid properties
rho_fluid = 1000
viscostity_din = 0.8e-03
viscostity_kin = viscostity_din / rho_fluid
termCondutivity_fluid = 0.6089
spHeat_fluid = 4137.9
termDiffusivity_fluid = termCondutivity_fluid / (rho_fluid * spHeat_fluid)
termDiffusivity_fluid = 1

# -------------------------------------------------------
#     Initial Condition
# -------------------------------------------------------

# Flow
vz = sp.zeros(nodes_global)
vr = sp.zeros(nodes_global)

# Heat
temp_old = sp.zeros(nodes_global)
q = 1
Q = sp.ones(nodes_global) * q

# -------------------------------------------------------
#     Matrix Assembly
# -------------------------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):

        v = [_ien[elem, 0], _ien[elem, 1], _ien[elem, 2]]
        axisym_tri.getAxiSym(v)

        ele_radius = (_y[v[0]] + _y[v[1]] + _y[v[2]])/3.

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]

                k_global[i_global, j_global]      += ele_radius*axisym_tri.kxx[i_local, j_local]+\
                                                     ele_radius*axisym_tri.kyy[i_local, j_local]
                mr_global[i_global, j_global]     += ele_radius*axisym_tri.mass[i_local, j_local]


    return k_global, mr_global


K, M = fem_matrix(x_global, y_global, num_ele_global, nodes_global, ien_global)

Mdt = M / dt

# ---------------------------------------
# Boundary and Initial Condition
# ---------------------------------------

dirichlet_len = len(global_mesh.dirichlet_points)


temp_a = sp.zeros(nodes_global)
dpdx = -2
for i in range(nodes_global):
    temp_a[i] = 1


# --------------------------------------
# LHS matrix with Dirichlet BC


LHS = Mdt + K * termDiffusivity_fluid
LHS_og = sp.copy(LHS)

cct = sp.zeros(nodes_global)
for i in range(dirichlet_len):
    index = int(global_mesh.dirichlet_points[i][0] - 1)
    value = global_mesh.dirichlet_points[i][1]
    temp_old[index] = value
    for j in range(nodes_global):
        cct[j] -= value * LHS_og[j, index]
        if j != index:
            LHS[index, j] = 0
            LHS[j, index] = 0
        else:
            LHS[index, j] = 1


temp_ana = sp.zeros(nodes_global)

from scipy import special

# maxsum = 1000
# for i in range(nodes_global):
#     z = x_global[i]
#     r = y_global[i]
#     for n in range(1,maxsum):
#         ln = n*(sp.pi/5.)
#         i0 = sp.special.iv(0, ln)
#         ir = sp.special.iv(0, r*ln)
#         c1 = (2./ln) * (sp.sin(2.5*ln))**2
#         c2 = 2.5 - sp.sin(10*ln)/(4*ln)
#         cn = c1 / (i0 * c2)
#         temp_ana[i] += cn * ir * sp.sin(z*ln)

def gz(_z):
    return (1 - (q*L**2)*(1./(2.*k)) * ((_z/L) - (_z/L)**2)) * sp.sin(ln * _z)

from scipy.integrate import quad
L = max(x_global)
k = 1
maxsum = 100
for i in range(nodes_global):
    print i, " / ", nodes_global
    z = x_global[i]
    r = y_global[i]
    for n in range(1, maxsum):
        ln = n*(sp.pi/5.)
        i0 = sp.special.iv(0, ln)
        ir = sp.special.iv(0, r*ln)
        c1 = quad(gz, 0, L)[0]
        c2 = 2.5 - sp.sin(10*ln)/(4*ln)
        cn = c1 / (i0 * c2)
        temp_ana[i] += cn * ir * sp.sin(z*ln)
    temp_ana[i] += (q*L**2)*(1./(2.*k)) * ((z/L) - (z/L)**2)


# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

for t in range(0, time):
    print "Solving System " + str((float(t)/time)*100) + "%"

    RHS = sp.dot(Mdt, temp_old) + cct + sp.dot(M, Q)
    for i in range(dirichlet_len):
        index = int(global_mesh.dirichlet_points[i][0]-1)
        RHS[index] = global_mesh.dirichlet_points[i][1]

    temp = sp.linalg.solve(LHS, RHS)

    error = temp_ana - temp
    error = sp.sqrt(sp.dot(error, error))
    print error


    # Salvar VTK
    vtk_t = IO.InOut(x_global, y_global, ien_global, nodes_global, num_ele_global, temp_old, None, None
                     , temp_ana, abs(temp_ana-temp), None, None)
    vtk_t.saveVTK(cwd+"/results", savename + str(t+1))

    if max(abs(temp - temp_a)) <= 1e-06:
        break

    temp_old = sp.copy(temp)

# ----------------- Fim de Loop -------------------
# -------------------------------------------------
