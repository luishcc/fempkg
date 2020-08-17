# -*- coding: utf-8 -*-

import os
import scipy as sp
from scipy import linalg
import InOut as Io
import meshio
import semiLagrangean as sl
from functions import *

#import GMesh as Gm                                      # No Cython
from cython_func.assembly import GMesh as Gm             # With Cython
from cython_func.assembly.cyassembly import fem_matrix   # With Cython

cwd = os.getcwd()

arquivo = "vivC"
sim_case = 'vivMovingTest'

import importlib.util
spec = importlib.util.spec_from_file_location("resultsdir",
        "/home/luis/fempkg/fem-1.0/resultsdir.py")
dir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dir)
results_path = dir.make_dir(sim_case)
vtk_path = results_path + '/vtk'


malha = Gm.GMesh("mesh/"+arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN
nodes = len(x)
num_ele = len(ien)

dt = 0.001
tempo = 400

p_smooth = 0.6
p_exp = 1.3

# ---------------------------------------
# Cylinder movement function
# ---------------------------------------

def  move_cylinder2(_nodes, _y, _y_max, _f_0, _t, _dt):
  vel = 2 * sp.pi * _f_0 * _y_max * sp.cos(2*sp.pi*_f_0*_t)
  _center = 0
  _len_cylinder = len(_nodes)
  for i in _nodes:
    _y[i] = _y[i] + vel*_dt
    _center += _y[i] / _len_cylinder
  return _center, vel , _y

# ---------------------------------------
# Wz, Psi and Initial Velocity
# ---------------------------------------

Wz = sp.zeros(nodes, dtype="float64")

# ---------------------------------------
# Matrices Assembly
# ---------------------------------------

Boundary = sp.zeros(nodes)
lista = []

cylinder = Gm.get_boundary_with_tag("mesh/"+arquivo+".msh", 5)

center = 0
len_cylinder = len(cylinder)

for i in cylinder:
    center += y[i] / len_cylinder

for i in range(nodes):
    if y[i] == 0.0 :
        Boundary[i] = i
    elif y[i] == 10.0 :
        Boundary[i] = i
    elif i in cylinder:
        Boundary[i] = i
    elif x[i] == 0.0:
        Boundary[i] = i
    elif x[i] == 32.5:
        Boundary[i] = i
    else:
        lista.append(i)

Boundary = sp.delete(Boundary, lista, axis=0)
num_bc = len(Boundary)


# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

#neighbour_ele, neighbour_nodes = sl.neighbourElements2(nodes, ien)


f = open('vivC-neighbourNodes.txt', 'r')
neighbour_nodes = []
for line in f:
    temp = line[1:-2].split(', ')
    for i in range(len(temp)):
        temp[i] = int(temp[i])
    neighbour_nodes.append(temp)
f.close()

f = open('vivC-neighbourElements.txt', 'r')
neighbour_ele = []
for line in f:
    temp = line[1:-2].split(', ')
    for i in range(len(temp)):
        temp[i] = int(temp[i])
    neighbour_ele.append(temp)
f.close()


# def get_distance_y(_y, _center):
#     size = len(_y)
#     distance = sp.zeros(size)
#     for i in range(size):
#         distance[i] = abs(_y[i] - _center)
#     return distance

def set_mesh_velocity(_y, _vel_c, _center):
    size = len(_y)
    _vv = sp.zeros(size)
    for i in  range(size):
        d = -1*abs(_y[i]-_center)
        _vv[i] = (_vel_c +_vel_c*sp.exp(-5))*sp.exp(d) - _vel_c*sp.exp(-5)
    return _vv

vx_mesh = sp.zeros(nodes)

for t in range(0, tempo-1):

    print()
    print("Solving System " + str((float(t)/(tempo-1))*100) + "%")

    cyl_vel = 0
    cylinder_center, cyl_vel, y = move_cylinder2(cylinder, y, 0.3, 16, t*dt, dt)

    vy_exp = set_mesh_velocity(y, cyl_vel, cylinder_center)

    vx_smooth, vy_smooth = Gm.weighted_smoothMesh(neighbour_nodes, Boundary,
                                                    x, y, dt)

    vx_mesh = p_smooth * vx_smooth
    vy_mesh = p_smooth * vy_smooth + p_exp * vy_exp

    for i in range(num_bc):
        index = int(Boundary[i])
        vx_mesh[index] = 0
        vy_mesh[index] = 0

        # if index in cylinder:
        #     vy_mesh[index] = cyl_vel

    x = x + vx_mesh * dt
    y = y + vy_mesh * dt


    # Save VTK
    vtk = Io.InOut(x, y, ien, len(x), len(ien), Wz, Wz, Wz,
                    None, None, Wz, Wz, vx_mesh, vy_exp)
    vtk.saveVTK(vtk_path, arquivo + '-' + str(t+1))
    vtk.saveVTK(vtk_path, arquivo + '-last')
