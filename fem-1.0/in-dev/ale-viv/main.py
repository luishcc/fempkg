# -*- coding: utf-8 -*-

import os

#os.environ["MKL_NUM_THREADS"] = "8"
#os.environ["NUMEXPR_NUM_THREADS"] = "8"
#os.environ["OMP_NUM_THREADS"] = "8"

import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import meshio
import semiLagrangean as sl
from timeit import default_timer as timer
from functions import *


time_start = timer()
cwd = os.getcwd()

arquivo = "vivTest"
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

dt = 0.001
tempo = 400
Re = 100
v_in = Re
psi_top = v_in * 10

p_lagrange = 0.0
p_smooth = 0.7
p_wave = 0.0

param = {'Re':Re, 'dt':dt, 'tempo':tempo, 'p_lagrange':p_lagrange, \
            'p_smooth':p_smooth, 'Mesh File':"mesh/"+arquivo+".msh"}
dir.sim_info_files(results_path, param)


exit()

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

# ---------------------------------------
# Wz, Psi e velocidade inicial
# ---------------------------------------

Psi_new = sp.zeros(nodes, dtype="float64")
Wz_new = sp.zeros(nodes, dtype="float64")
vx = sp.zeros(nodes, dtype="float64") + v_in
vy = sp.zeros(nodes, dtype="float64")

# ---------------------------------------
# Montagem de matrizes
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

# ---------------------------------------
# Condições de contorno e Inicial
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
Wz_old = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx))) * w_temp

K_psi = sp.copy(K)
ccpsi = sp.zeros(nodes)

for i in range(num_bc):
    index = int(Boundary[i])
    value = psi_bc[index]
    if y[index] == 0.0 or y[index] == 10.0 or \
    (index in cylinder) or x[index] == 0:
        ccpsi[index] -= value * K[j, index]
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
neighbour_ele, neighbour_nodes = sl.neighbourElements2(nodes, ien)
time_end_neighbour = timer()
print("Neighbour structures: ", time_end_neighbour - time_start_neighbour)


time_avg_loop = 0
for t in range(0, tempo-1):

    time_start_loop = timer()
    print()
    print("Solving System " + str((float(t)/(tempo-1))*100) + "%")

    time_start_smooth = timer()

    cyl_vel = 0
    y, cylinder_center, cyl_vel = move_cylinder(cylinder, y, 0.3, 16, t*dt, dt)
    for i in cylinder:
        psi_bc[i] = 0.1 * psi_top * cylinder_center

#    vx_smooth, vy_smooth = Gm.smoothMesh(neighbour_nodes, malha, x, y, dt)
    vx_smooth, vy_smooth = Gm.weighted_smoothMesh(neighbour_nodes, Boundary,
                                                    x, y, dt)

    time_end_smooth = timer()
    print("Smoothing time: ", time_end_smooth - time_start_smooth)


#    vx_wave = sp.random.rand(nodes) * 0.2*sp.cos(2*sp.pi*5*t*dt)
#    vy_wave = sp.random.rand(nodes) * 0.2*sp.cos(2*sp.pi*5*t*dt)


    time_start_ale = timer()
    vxAle = p_lagrange * vx + p_smooth * vx_smooth #+ p_wave * vx_wave
    vyAle = p_lagrange * vy + p_smooth * vy_smooth #+ p_wave * vy_wave

    for i in range(num_bc):
        index = int(Boundary[i])
        vxAle[index] = 0
        vyAle[index] = 0

        if x[index] == 0 or y[index] == 10 or y[index] == 0:
            vx[index] = v_in
            vy[index] = 0
        if index in cylinder:
            vx[index] = 0
            #vy[index] = 0
            vy[index] = cyl_vel #(Change to displacement velocity?)


    vx_sl = vx - vxAle
    vy_sl = vy - vyAle
    Wz_dep = sl.Linear2D(nodes, neighbour_ele, ien, x, y, vx_sl,
                            vy_sl, dt, Wz_old)

    # Solução de Wz e Psi
    x = x + vxAle * dt
    y = y + vyAle * dt
    K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

    Minv = linalg.inv(M)

    time_end_ale = timer()
    print("ALE step time: ", time_end_ale - time_start_ale)


    time_start_wz = timer()

    # B.C. Vorticidade
    Wcc = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx))) * w_temp
    ccomega = sp.zeros(nodes)


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

    time_end_wz = timer()
    print("Omega solution time: ", time_end_wz - time_start_wz)

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

    # Salvar VTK
    vtk = Io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, Wz_dep,
                    None, None, vx, vy, vx_sl, vy_sl)
    vtk.saveVTK(vtk_path, arquivo + '-' + str(t+1))

    Psi_old = sp.copy(Psi_new)
    Wz_old = sp.copy(Wz_new)

    # Calculo de Vx e Vy
    vx = sp.dot(Minv, sp.dot(Gy, Psi_new))
    vy = -1.0 * sp.dot(Minv, sp.dot(Gx, Psi_new))

    time_end_loop = timer()
    time_avg_loop += time_end_loop - time_start_loop
    print("Iteration time: ", time_end_loop - time_start_loop)
    print("Estimated time to finish simulation: ",
            (time_end_loop - time_start_loop)* (tempo-t+1) )

#----------------- Fim de Loop -------------------
#-------------------------------------------------
time_end = timer()
print()
print("Iteration average time: ", time_avg_loop/(t+1))
print("Total runtime: ", time_end - time_start)
