# -*- coding: utf-8 -*-

import os

os.environ["MKL_NUM_THREADS"] = "5"
os.environ["NUMEXPR_NUM_THREADS"] = "5"
os.environ["OMP_NUM_THREADS"] = "5"

import scipy as sp
import numpy as np
from scipy import linalg

import lib.inOut as io
import lib.semiLagrangian as sl
import lib.mesh as m
from lib.cython_func.assembly.assembly import fem_matrix   # With Cython
# from lib.assembly import fem_matrix   # No Cython
from lib.apply_bc import *
import lib.resultsdir as rt_save

from timeit import default_timer as timer


time_start = timer()
cwd = os.getcwd()

start_file = 'last-200fine01'
start_from_file = True
#start_from_file = False

msh_file = "cylinder-fine"
sim_case = 'flowAroundCylinder'
sim_type='fixed'
#sim_type='moving'

results_path = rt_save.make_dir(sim_case)
vtk_path = results_path + '/vtk'

print('-------------------------------------------------------------')
print()
print('  Starting Simulation: {}'.format(sim_case))
print('  Saving Results in: {}'.format(results_path))
print()
print('-------------------------------------------------------------')

time_start_read = timer()

mesh = m.Mesh("mesh/{}/{}.msh".format(sim_case,msh_file))
x = mesh.x
y = mesh.y
ien = mesh.ien
NN = mesh.num_nodes
NE = mesh.num_elem

time_end_read = timer()

print("Read .msh in: ", time_end_read - time_start_read)

dt = 0.008
steps = 10000
vtk_steps = 1

Re = 200
v_in = 1
psi_top = max(y)

p_lagrange = 0.0
p_smooth = 0.0
p_exp = 1.2

param = {'Re':Re, 'dt':dt, 'Steps':steps, 'p_lagrange':p_lagrange, \
            'p_smooth':p_smooth, 'p_exp':p_exp ,\
            'Mesh File':"mesh/{}/{}.msh".format(sim_case,msh_file)}
rt_save.sim_info_files(results_path, param)


# ---------------------------------------
# Mesh and Cylinder movement functions
# ---------------------------------------

def set_mesh_velocity(_y, _vel_c, _center):
    size = len(_y)
    _vv = sp.zeros(size)
    for i in  range(size):
        d = -1*abs(_y[i]-_center)
        _vv[i] = (_vel_c +_vel_c*sp.exp(-5))*sp.exp(d) - _vel_c*sp.exp(-5)
    return _vv

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


if start_from_file:
    import meshio
    temp_file = meshio.read('last/'+start_file+'.vtk')
    vx = temp_file.point_data['Velocity'][:,0]
    vy = temp_file.point_data['Velocity'][:,1]
    Psi_new = temp_file.point_data['Psi'][:,0]
    Wz_old = temp_file.point_data['Omega'][:,0]



# ---------------------------------------
# Boundary and Initial Conditions
# ---------------------------------------

time_start_bc = timer()

psi_dirichlet_nodes = mesh.set_dirichlet_nodes(['top', 'bot',
                                                'inlet', 'cylinder'])
psi_bc_value = np.zeros(NN)
cylinder_nodes = mesh.get_boundary_with_name('cylinder')

omega_null_bc = np.ones(NN)
omega_bc_value = np.zeros(NN)

for i in mesh.boundary_nodes:
    vy[i] = 0

for i in mesh.get_boundary_with_name('top'):
    psi_bc_value[i] = psi_top
    vx[i] = v_in
    omega_null_bc[i] = 0

for i in mesh.get_boundary_with_name('bot'):
    psi_bc_value[i] = 0
    vx[i] = v_in
    omega_null_bc[i] = 0

for i in cylinder_nodes:
    psi_bc_value[i] = 0.5 * psi_top
    vx[i] = 0
    vy[i] = 0

for i in mesh.get_boundary_with_name('inlet'):
    psi_bc_value[i] = 0.1 * y[i] * psi_top
    vx[i] = v_in
    omega_null_bc[i] = 0

num_bc = len(mesh.boundary_nodes)


time_end_bc = timer()
print("Boundary and initial conitions: ", time_end_bc - time_start_bc)

# ----------------------------------------------------------
# ---------------------- Time Loop -------------------------

# exit()

v_c = np.zeros(NN)
time_avg_loop = 0
iter = 10000

K, M, Gx, Gy = fem_matrix(mesh)
Minv = linalg.inv(M)

for t in range(0, int(steps/vtk_steps)):
    noise = sp.random.rand(NN)* v_in * 0.0005
    for tt in range(vtk_steps):
        iter += 1
        if iter <=20:
            vy = vy * noise
            vx = vx * noise
        time_start_loop = timer()
        print()
        print("Solving System " + str((float(iter/10000)/(steps-1))*100) + "%")

        cyl_vel = 0


        Wz_dep = sl.Linear2D(NN, mesh.neighbour_elements, ien, x, y, vx,
                                vy, dt, Wz_old)

        # Wz , Psi Solution

        time_start_wz = timer()

        # B.C. Vorticity
        omega_bc_value = np.dot(Minv, (np.dot(Gx, vy) - np.dot(Gy, vx))) * omega_null_bc


        # Solve vorticity Transport
        LHS = M / dt + K / Re
        LHS_omega, omega_bc_RHS = apply_bc_dirichlet(LHS, mesh,
                                                     mesh.boundary_nodes,
                                                     omega_bc_value)

        F_omega = np.dot(M / dt, Wz_dep) + omega_bc_RHS /Re
        for i in mesh.boundary_nodes:
            F_omega[i] = omega_bc_value[i]

        Wz_new = sp.linalg.solve(LHS_omega, F_omega)

        time_end_wz = timer()
        print("Omega solution time: ", time_end_wz - time_start_wz)

        # Solve Stream Function
        time_start_psi = timer()
        K_psi, psi_bc_RHS = apply_bc_dirichlet(K, mesh, psi_dirichlet_nodes,
                                                psi_bc_value)
        F_psi = np.dot(M, Wz_new) + psi_bc_RHS
        for i in psi_dirichlet_nodes:
            F_psi[i] = psi_bc_value[i]
        Psi_new = sp.linalg.solve(K_psi, F_psi)
        time_end_psi = timer()
        print("Psi solution time: ", time_end_psi - time_start_psi)

        # Saving last step solution
        Psi_old = np.copy(Psi_new)
        Wz_old = np.copy(Wz_new)

        for i in mesh.boundary_nodes:
            if i in mesh.get_boundary_with_name('inlet'):
                vx[i] = v_in
                vy[i] = 0
            if i in mesh.get_boundary_with_name('top') or \
               i in mesh.get_boundary_with_name('bot'):
                vx[i] = v_in
                vy[i] = 0
            if i in cylinder_nodes:
                vx[i] = 0
                vy[i] = 0

        LHS_vx, BC_vx = apply_bc_dirichlet(M, mesh, psi_dirichlet_nodes, vx)
        LHS_vy, BC_vy = apply_bc_dirichlet(M, mesh, psi_dirichlet_nodes, vy)
        F_vy = np.dot(Gx, Psi_new) + BC_vy
        F_vx = np.dot(Gy, Psi_new) + BC_vx
        for i in psi_dirichlet_nodes:
            F_vx[i] = vx[i]
            F_vy[i] = vy[i]
        vx = sp.linalg.solve(LHS_vx, F_vx)
        vy = -1*sp.linalg.solve(LHS_vy, F_vy)


        time_end_loop = timer()
        time_avg_loop += time_end_loop - time_start_loop
        print("Iteration time: ", time_end_loop - time_start_loop)
        print("Estimated time to finish simulation: ",
                (time_end_loop - time_start_loop)* (steps-t+1) )

    # Save VTK
    vtk = io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, Wz_dep,
                    None, None, vx, vy, vx, vy)
    vtk.saveVTK(vtk_path,  'vtk'+str(iter))
    vtk.saveVTK(vtk_path, 'last')

#-----------------  Loop End  --------------------
#-------------------------------------------------
time_end = timer()
print()
print("Iteration average time: ", time_avg_loop/(t+1))
print("Total runtime: ", time_end - time_start)
