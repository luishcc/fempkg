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

print "neighbour"
neighbour = sl.neighbourElements(nodes, ien)

for t in range(iterations):

    print t, " / ", time
    zk_new = sl.Linear2D(nodes, neighbour, ien, x, y, vx, vy, dt, zk)

    # Salvar VTK
    vtk_t = IO.InOut(x, y, ien, nodes, num_ele, zk, None, None
                     , None, None, vx, vy)
    vtk_t.saveVTK(cwd+"/results", mesh_file + str(t))

    zk = sp.copy(zk_new)

# ----------------- Fim de Loop -------------------
# -------------------------------------------------
