# -*- coding: utf-8 -*-

import os

os.environ["MKL_NUM_THREADS"] = "3"
os.environ["NUMEXPR_NUM_THREADS"] = "3"
os.environ["OMP_NUM_THREADS"] = "3"

import scipy as sp
from scipy import linalg
import InOut as Io
import meshio
import semiLagrangean as sl
from timeit import default_timer as timer
from functions import *

#import GMesh as Gm                                      # No Cython
from cython_func.assembly import GMesh as Gm             # With Cython
from cython_func.assembly.cyassembly import fem_matrix   # With Cython



time_start = timer()
cwd = os.getcwd()

arquivo = "vivC"
sim_case = 'vivMovingCylinder'


import importlib.util
spec = importlib.util.spec_from_file_location("resultsdir",
        "/home/luis/fempkg/fem-1.0/resultsdir.py")
dir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dir)
results_path = dir.make_dir(sim_case)
vtk_path = results_path + '/vtk'

print('-------------------------------------------------------------')
print()
print('  Starting Simulation: {}'.format(sim_case))
print('  Saving Results in: {}'.format(results_path))
print()
print('-------------------------------------------------------------')

time_start_read = timer()

malha = Gm.GMesh("mesh/"+arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN
nodes = len(x)
num_ele = len(ien)

time_end_read = timer()

print("Read .msh in: ", time_end_read - time_start_read)

dt = 0.05
tempo = 2000
Re = 10
v_in = 1
psi_top = 10

p_lagrange = 0.0
p_smooth = 0.0
p_exp = 1.2

param = {'Re':Re, 'dt':dt, 'tempo':tempo, 'p_lagrange':p_lagrange, \
            'p_smooth':p_smooth, 'p_exp':p_exp ,\
            'Mesh File':"mesh/"+arquivo+".msh"}
dir.sim_info_files(results_path, param)



# ---------------------------------------
# Cylinder movement function
# ---------------------------------------

def  move_cylinder(_nodes, _y, _y_max, _f_0, _t, _dt):
  vel = 2 * sp.pi * _f_0 * _y_max * sp.cos(2*sp.pi*_f_0*_t)
  _center = 0
  _len_cylinder = len(_nodes)
  for i in _nodes:
    _y[i] = _y[i] + vel*_dt
    _center += _y[i] / _len_cylinder
  return _y, _center, vel

def  move_cylinder2(_nodes, _y, _y_max, _f_0, _t, _dt):
  vel = 2 * sp.pi * _f_0 * _y_max * sp.cos(2*sp.pi*_f_0*_t)
  _center = 0
  _len_cylinder = len(_nodes)
  for i in _nodes:
    _center += _y[i] / _len_cylinder
  return _center, vel

# ---------------------------------------
# Wz, Psi and Initial Velocity
# ---------------------------------------

Psi_new = sp.zeros(nodes, dtype="float64")
Wz_new = sp.zeros(nodes, dtype="float64")
vx = sp.zeros(nodes, dtype="float64") #+ v_in
vy = sp.zeros(nodes, dtype="float64")

# ---------------------------------------
# Matrices Assembly
# ---------------------------------------

time_start_assembly = timer()

K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

time_end_assembly = timer()
print("First assemble in: ", time_end_assembly - time_start_assembly)


time_start_bc = timer()
Boundary = sp.zeros(nodes)
lista = []

psi_bc = sp.zeros(nodes)

cylinder = Gm.get_boundary_with_tag("mesh/"+arquivo+".msh", 5)
center = 0
len_cylinder = len(cylinder)

for i in cylinder:
    center += y[i] / len_cylinder

for i in range(nodes):
    if y[i] == 0.0 :
        Boundary[i] = i
        psi_bc[i] = 0
    elif y[i] == 10.0 :
        Boundary[i] = i
        psi_bc[i] = psi_top
    elif i in cylinder:
        Boundary[i] = i
        psi_bc[i] = 0.5 * psi_top
    elif x[i] == 0.0:
        psi_bc[i] = 0.1 * psi_top * y[i]
        Boundary[i] = i
    elif x[i] == 32.5:
        Boundary[i] = i
    else:
        lista.append(i)

Boundary = sp.delete(Boundary, lista, axis=0)
num_bc = len(Boundary)

# cylinder = sp.ones(num_bc)*-1

# ---------------------------------------
# Boundary and Initial Conditions
# ---------------------------------------

w_temp = sp.ones(nodes)
bc_omega = sp.zeros(nodes)
for i in Boundary:
    j = int(i)
    if y[j] == 10 or y[j] == 0 or x[j] == 0:
        vx[j] = v_in
        w_temp[j] = 0
        vy[j] = 0.
    if j in cylinder:
        vx[j] = 0
        vy[j] = 0


Minv = sp.linalg.inv(M)
#Wz_old = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx))) * w_temp
Wz_old = sp.zeros(nodes)

K_psi = sp.copy(K)
ccpsi = sp.zeros(nodes)

for i in range(num_bc):
    index = int(Boundary[i])
    value = psi_bc[index]
    if y[index] == 0.0 or y[index] == 10.0 or \
    (index in cylinder) or x[index] == 0:
        ccpsi[index] -= value * K[i, index]
        for j in range(nodes):
            if j != index:
                K_psi[index, j] = 0
                K_psi[j, index] = 0
            else:
                K_psi[index, j] = 1

F_psi = sp.dot(M, Wz_old) + ccpsi
for i in range(num_bc):
    index = int(Boundary[i])
    F_psi[index] = psi_bc[index]

Psi_old = sp.linalg.solve(K_psi, F_psi)

time_end_bc = timer()

print("Boundary and initial conitions: ", time_end_bc - time_start_bc)
sp.random.seed(1)

# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

time_start_neighbour = timer()
# neighbour_ele, neighbour_nodes = sl.neighbourElements2(nodes, ien)


f = open(arquivo+'-neighbourNodes.txt', 'r')
neighbour_nodes = []
for line in f:
    temp = line[1:-2].split(', ')
    for i in range(len(temp)):
        temp[i] = int(temp[i])
    neighbour_nodes.append(temp)
f.close()

f = open(arquivo+'-neighbourElements.txt', 'r')
neighbour_ele = []
for line in f:
    temp = line[1:-2].split(', ')
    for i in range(len(temp)):
        temp[i] = int(temp[i])
    neighbour_ele.append(temp)
f.close()

time_end_neighbour = timer()
print("Neighbour structures: ", time_end_neighbour - time_start_neighbour)

def set_mesh_velocity(_y, _vel_c, _center):
    size = len(_y)
    _vv = sp.zeros(size)
    for i in  range(size):
        d = -1*abs(_y[i]-_center)
        _vv[i] = (_vel_c +_vel_c*sp.exp(-5))*sp.exp(d) - _vel_c*sp.exp(-5)
    return _vv

v_c = sp.zeros(nodes)
time_avg_loop = 0
for t in range(0, tempo-1):

    time_start_loop = timer()
    print()
    print("Solving System " + str((float(t)/(tempo-1))*100) + "%")

    time_start_smooth = timer()

    cyl_vel = 0
#    y, cylinder_center, cyl_vel = move_cylinder(cylinder, y, 0.3, 16, t*dt, dt)
    # cylinder_center, cyl_vel = move_cylinder2(cylinder, y, 0.3, 16, t*dt, dt)
#    for i in cylinder:
#        psi_bc[i] = 0.1 * psi_top * cylinder_center

#    vy_exp = set_mesh_velocity(y, cyl_vel, cylinder_center)

#    vx_smooth, vy_smooth = Gm.smoothMesh(neighbour_nodes, malha, x, y, dt)
    vx_smooth, vy_smooth = Gm.weighted_smoothMesh(neighbour_nodes, Boundary,
                                                    x, y, dt)

    time_end_smooth = timer()
    print("Smoothing time: ", time_end_smooth - time_start_smooth)


    time_start_ale = timer()
    vx_mesh = p_lagrange * vx + p_smooth * vx_smooth
    vy_mesh = p_lagrange * vy + p_smooth * vy_smooth #+ p_exp * vy_exp

    for i in range(num_bc):
        index = int(Boundary[i])
        vx_mesh[index] = 0
        vy_mesh[index] = 0

        # if x[index] == 0 or y[index] == 10 or y[index] == 0:
        #     vx[index] = v_in
        #     vy[index] = 0
        # if index in cylinder:
        #     vx[index] = 0
        #     vy[index] = 0
#            vy[index] = cyl_vel
#            vy_mesh[index] = cyl_vel
#            v_c[index] = cyl_vel


    x = x + vx_mesh * dt
#    y = y + (vy_mesh - v_c) * dt
    y = y + vy_mesh  * dt

    vxAle = vx - vx_mesh
    vyAle = vy - vy_mesh

    Wz_dep = sl.Linear2D(nodes, neighbour_ele, ien, x, y, vxAle,
                            vyAle, dt, Wz_old)

    # Wz , Psi Solution

    K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

    #Minv = linalg.inv(M)

    time_end_ale = timer()
    print("ALE step time: ", time_end_ale - time_start_ale)


    time_start_wz = timer()

    # B.C. Vorticity
    Wcc = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx))) * w_temp
    ccomega = sp.zeros(nodes)

    # Solve vorticity Transport
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

    F_omega = sp.dot(M / dt, Wz_dep) + ccomega # / Re

    for i in range(num_bc):
        index = int(Boundary[i])
        F_omega[index] = Wcc[index]

    Wz_new = sp.linalg.solve(LHS_omega, F_omega)

    time_end_wz = timer()
    print("Omega solution time: ", time_end_wz - time_start_wz)

    # Solve Stream Function
    time_start_psi = timer()

    K_psi = sp.copy(K)
    ccpsi = sp.zeros(nodes)

    for i in range(num_bc):
        index = int(Boundary[i])
        value = psi_bc[index]
        if y[index] == 0.0 or y[index] == 10.0 or \
        (index in cylinder) or x[index] == 0:
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
        F_psi[index] = psi_bc[index]

    Psi_new = sp.linalg.solve(K_psi, F_psi)

    time_end_psi = timer()
    print("Psi solution time: ", time_end_psi - time_start_psi)

    # Check for steady state with vorticity solution
    err = Wz_new - Wz_old
    err = sp.sqrt(sp.dot(err,err))
    if err < 1e-8:
        break

    # Saving last step solution
    Psi_old = sp.copy(Psi_new)
    Wz_old = sp.copy(Wz_new)

    # Calculate Vx e Vy
#    vx = sp.dot(Minv, sp.dot(Gy, Psi_new))
#    vy = -1.0 * sp.dot(Minv, sp.dot(Gx, Psi_new))
    vx = sp.linalg.solve(M, sp.dot(Gy, Psi_new))
    vy = -1.0 * sp.linalg.solve(M, sp.dot(Gx, Psi_new))
    # Setting velocity BC
    for i in range(num_bc):
        index = int(Boundary[i])
        if x[index] == 0 :
            vx[index] = v_in
            vy[index] = 0
        if y[index] == 10 or y[index] == 0:
            vx[index] = v_in
            vy[index] = 0
        if index in cylinder:
            vx[index] = 0
            vy[index] = cyl_vel


    # Save VTK
    vtk = Io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, Wz_dep,
                    None, None, vx, vy, vxAle, vyAle)
    vtk.saveVTK(vtk_path, arquivo + '-' + str(t+1))
    vtk.saveVTK(vtk_path, arquivo + '-last')

    time_end_loop = timer()
    time_avg_loop += time_end_loop - time_start_loop
    print("Iteration time: ", time_end_loop - time_start_loop)
    print("Estimated time to finish simulation: ",
            (time_end_loop - time_start_loop)* (tempo-t+1) )

#-----------------  Loop End  --------------------
#-------------------------------------------------
time_end = timer()
print()
print("Iteration average time: ", time_avg_loop/(t+1))
print("Total runtime: ", time_end - time_start)
