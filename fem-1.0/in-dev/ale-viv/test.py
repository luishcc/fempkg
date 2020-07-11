# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import meshio
import os
import semiLagrangean as sl

cwd = os.getcwd()

arquivo = "vivC"

malha = Gm.GMesh("mesh/"+arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN
nodes = len(x)
num_ele = len(ien)

dt = 0.0001
tempo = 300

p_smooth = 0.7



# ---------------------------------------
# Forced Cylinder Vibration
# ---------------------------------------

def move_cylinder(_nodes, _y, _y_max, _f_0, _t, _dt):
    vel = 2*sp.pi*_f_0*_y_max*sp.cos(2*sp.pi*_f_0*_t)
    for i in _nodes:
        _y[i] = _y[i] + vel*_dt
    return _y


# ---------------------------------------
# Wz, Psi e velocidade inicial
# ---------------------------------------

Psi_old = sp.zeros(nodes)


# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

Boundary = sp.zeros(nodes)
lista = []

cylinder = Gm.get_boundary_with_tag("mesh/"+arquivo+".msh", 5)


for i in range(nodes):
    if y[i] == 0.0 :
        Boundary[i] = i
    elif y[i] == 10.0 :
        Boundary[i] = i
    elif i in cylinder:
        Boundary[i] = i
    elif x[i] == 0.0 or x[i] == 32.5:
        Boundary[i] = i
    else:
        lista.append(i)

Boundary = sp.delete(Boundary, lista, axis=0)
num_bc = len(Boundary)

sp.random.seed(1)

# -------------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

print("Neighbour Elements ")
neighbour_ele, neighbour_nodes = sl.neighbourElements2(nodes, ien)

neighbourfile = open(arquivo + "-neighbourNodes.txt", "w")
neighbourfile2 = open(arquivo + "-neighbourElements.txt", "w")
for i in range(len(neighbour_nodes)):
    nodes_temp = str(neighbour_nodes[i])
    ele_temp = str(neighbour_ele[i])
    neighbourfile.writelines(nodes_temp+"\n")
    neighbourfile2.writelines(ele_temp+"\n")

# for i in range(1):
#     vx_smooth, vy_smooth = Gm.weighted_smoothMesh(neighbour_nodes, Boundary, x, y, dt)
#     x = x + p_smooth*vx_smooth  * dt
#     y = y + p_smooth*vy_smooth * dt

for t in range(0, tempo-1):
    print("Solving System " + str((float(t)/(tempo-1))*100) + "%")

    y = move_cylinder(cylinder, y, 0.3, 16, t*dt, dt)


    vx_smooth, vy_smooth = Gm.weighted_smoothMesh(neighbour_nodes, Boundary, x, y, dt)


    vxAle = p_smooth * vx_smooth
    vyAle = p_smooth * vy_smooth


    for i in range(num_bc):
        index = int(Boundary[i])
        vxAle[index] = 0
        vyAle[index] = 0


    x = x + vxAle  * dt
    y = y + vyAle * dt

    # Salvar VTK
    vtk = Io.InOut(x, y, ien, len(x), len(ien), Psi_old, Psi_old, Psi_old,
                    None, None, Psi_old, Psi_old, Psi_old, Psi_old)
    vtk.saveVTK(cwd+"/results", "move" + str(t+1))



#----------------- Fim de Loop -------------------
#-------------------------------------------------

# for j in range(1000):
#     vx_temp, vy_temp, x, y = Gm.smoothMesh_v2(neighbour_nodes, Boundary,
#                                               x, y, dt)
#
# vtk = Io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, Wz_dep,
#                 None, None, vx, vy, vx_sl, vy_sl)
# vtk.saveVTK(cwd+"/results", "move" + str(t+2))
