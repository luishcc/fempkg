# -*- coding: utf-8 -*-

import os

#os.environ["MKL_NUM_THREADS"] = "1"
#os.environ["NUMEXPR_NUM_THREADS"] = "1"
#os.environ["OMP_NUM_THREADS"] = "1"

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

start_file = 'last.vtk'
start_from_file = True
#start_from_file = False

msh_file = "vivC"
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

dt = 0.01
steps = 5000
vtk_steps = 1

Re = 120
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

Psi_new = np.zeros(NN, dtype="float64")
Wz_new = np.zeros(NN, dtype="float64")
vx = np.zeros(NN, dtype="float64")
vy = np.zeros(NN, dtype="float64")

if start_from_file:
    import meshio
    temp_file = meshio.read(start_file)
    vx = temp_file.point_data['Velocity'][:,0]
    vy = temp_file.point_data['Velocity'][:,1]
    Psi_new = temp_file.point_data['Psi']
    vx = temp_file.point_data['Omega']

# ---------------------------------------
# Matrices Assembly
# ---------------------------------------

time_start_assembly = timer()

# K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien) # No Cython
K, M, Gx, Gy = fem_matrix(mesh)

time_end_assembly = timer()
print("First assemble in: ", time_end_assembly - time_start_assembly)

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


Minv = sp.linalg.inv(M)
Wz_old = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx))) * omega_null_bc
#Wz_old = np.zeros(NN)


K_psi, psi_bc_RHS = apply_bc_dirichlet(K, mesh, psi_dirichlet_nodes,
                                        psi_bc_value)
F_psi = np.dot(M, Wz_old) + psi_bc_RHS
for i in psi_dirichlet_nodes:
    F_psi[i] = psi_bc_value[i]

Psi_old = sp.linalg.solve(K_psi, F_psi)

time_end_bc = timer()
print("Boundary and initial conitions: ", time_end_bc - time_start_bc)

# ----------------------------------------------------------
# ---------------------- Time Loop -------------------------

# exit()

v_c = np.zeros(NN)
time_avg_loop = 0
iter = 0
for t in range(0, int(steps/vtk_steps)):
    for tt in range(vtk_steps):
        iter += 1
        time_start_loop = timer()
        print()
        print("Solving System " + str((float(iter)/(steps-1))*100) + "%")

        cyl_vel = 0
        # if sim_type == 'moving':
        #     y, cylinder_center, cyl_vel = move_cylinder(cylinder,
        #                                                 y, 0.3, 16, t*dt, dt)
        #     # cylinder_center, cyl_vel = move_cylinder2(cylinder,
        #     #                                           y, 0.3, 16, t*dt, dt)
        #     for i in cylinder_nodes:
        #         psi_bc_value[i] = 0.1 * cylinder_center
        #     vy_exp = set_mesh_velocity(y, cyl_vel, cylinder_center)

        time_start_smooth = timer()
        vx_smooth, vy_smooth = m.weighted_smoothMesh(mesh.neighbour_nodes,
                                                     mesh.boundary_nodes,
                                                     x, y, dt)
        time_end_smooth = timer()
        print("Smoothing time: ", time_end_smooth - time_start_smooth)

        time_start_ale = timer()
        vx_mesh = p_lagrange * vx + p_smooth * vx_smooth
        vy_mesh = p_lagrange * vy + p_smooth * vy_smooth #+ p_exp * vy_exp

        for i in mesh.boundary_nodes:
            vx_mesh[i] = 0
            vy_mesh[i] = 0

#            if i in cylinder_nodes:
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

        Wz_dep = sl.Linear2D(NN, mesh.neighbour_elements, ien, x, y, vxAle,
                                vyAle, dt, Wz_old)

        # Wz , Psi Solution

        K, M, Gx, Gy = fem_matrix(mesh)

        Minv = linalg.inv(M)

        time_end_ale = timer()
        print("ALE step time: ", time_end_ale - time_start_ale)


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


        # Check for steady state with vorticity solution
        # err = Wz_new - Wz_old
        # err = np.sqrt(np.dot(err,err))
        # if err < 1e-8:
        #     exit()

        # Saving last step solution
        Psi_old = np.copy(Psi_new)
        Wz_old = np.copy(Wz_new)

        # Calculate Vx e Vy
    #    vx = sp.dot(Minv, sp.dot(Gy, Psi_new))
    #    vy = -1.0 * sp.dot(Minv, sp.dot(Gx, Psi_new))
        vx = sp.linalg.solve(M, np.dot(Gy, Psi_new))
        vy = -1.0 * sp.linalg.solve(M, np.dot(Gx, Psi_new))
        # Setting velocity BC
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
                vy[i] = cyl_vel

        time_end_loop = timer()
        time_avg_loop += time_end_loop - time_start_loop
        print("Iteration time: ", time_end_loop - time_start_loop)
        print("Estimated time to finish simulation: ",
                (time_end_loop - time_start_loop)* (steps-t+1) )

    # Save VTK
    vtk = io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, Wz_dep,
                    None, None, vx, vy, vxAle, vyAle)
    vtk.saveVTK(vtk_path,  'vtk'+str(iter))
    vtk.saveVTK(vtk_path, 'last')

#-----------------  Loop End  --------------------
#-------------------------------------------------
time_end = timer()
print()
print("Iteration average time: ", time_avg_loop/(t+1))
print("Total runtime: ", time_end - time_start)
