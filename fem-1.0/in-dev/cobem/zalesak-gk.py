# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
import semiLagrangean as sl
import os

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

mesh_file = "zalesak-coarse"

fluid_mesh = gm.GMesh("mesh/" + mesh_file + ".msh")

x = fluid_mesh.X
y = fluid_mesh.Y
ien = fluid_mesh.IEN
nodes = len(x)
num_ele = len(ien)


# -------------------------------------------------------
#     Simulation Parameters
# -------------------------------------------------------


omega = 0.5
iterations = 1000
time = (sp.pi*2) / omega
dt = (1./iterations) * time


ox = 2.0
oy = 2.0

# -------------------------------------------------------
#     Zalesak Conditions
# -------------------------------------------------------

R = 0.5
r = 0.4
s = 0.06
Cx = 2.0
Cy = 2.75

zk = sp.ones(nodes)
vx = sp.zeros(nodes)
vy = sp.zeros(nodes)
for i in range(nodes):
    vx[i] = -1 * omega * (y[i] - oy)
    vy[i] = omega * (x[i] - ox)

    if (Cx - 0.5*s < x[i] < Cx + 0.5*s) and (y[i] < Cy + (R - r)):
        zk[i] = 0
    if (x[i] - Cx)**2 + (y[i] - Cy)**2 > R**2:
        zk[i] = 0

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
        alfa = 0
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

Conv = sp.dot(sp.diag(vx), Gx) + sp.dot(sp.diag(vy), sp.diag(Gy))


LHS = M / dt + Conv

for t in range(iterations):
    print t, " / ", iterations

    RHS = sp.dot(M/dt, zk)
    zk_new = sp.linalg.solve(LHS, RHS)

    # Salvar VTK
    vtk_t = IO.InOut(x, y, ien, nodes, num_ele, zk, None, None
                     , None, None, vx, vy)
    vtk_t.saveVTK(cwd+"/resultsTest", mesh_file + "-GK" + str(t))

    zk = sp.copy(zk_new)

# ----------------- Fim de Loop -------------------
# -------------------------------------------------
