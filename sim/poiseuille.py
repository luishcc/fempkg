# -*- coding: utf-8 -*-

import os

#os.environ["MKL_NUM_THREADS"] = "3"
#os.environ["NUMEXPR_NUM_THREADS"] = "3"
#os.environ["OMP_NUM_THREADS"] = "3"

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

msh_file = "poiseuille-5"
sim_case = 'poiseuille'


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
steps = 1000
vtk_steps = 1

Re = 1
v_in = 1
psi_top = 1

p_lagrange = 0.0
p_smooth = 0.
p_w = 0.

param = {'Re':Re, 'dt':dt, 'Steps':steps, 'p_lagrange':p_lagrange, \
            'p_smooth':p_smooth, 'p_w':p_w ,\
            'Mesh File':"mesh/{}/{}.msh".format(sim_case,msh_file)}
rt_save.sim_info_files(results_path, param)


# ---------------------------------------
# Wz, Psi and Initial Velocity
# ---------------------------------------

Psi_new = np.zeros(NN, dtype="float64")
Wz_new = np.zeros(NN, dtype="float64")
vx = np.zeros(NN, dtype="float64")+1
vy = np.zeros(NN, dtype="float64")

vx_a = np.zeros(NN, dtype="float64")
w_a = np.zeros(NN, dtype="float64")
for i in range(NN):
    vx_a[i] = 6*y[i]*(1-y[i])
    w_a[i] = 12*y[i] - 6


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
                                                'inflow'])
psi_bc_value = np.zeros(NN)

omega_null_bc = np.ones(NN)
omega_bc_value = np.zeros(NN)

for i in mesh.boundary_nodes:
    vy[i] = 0

for i in mesh.get_boundary_with_name('top'):
    psi_bc_value[i] = psi_top
    vx[i] = 0
    omega_bc_value[i] = 6
    omega_null_bc[i] = 0


for i in mesh.get_boundary_with_name('bot'):
    psi_bc_value[i] = 0
    vx[i] = 0
    omega_bc_value[i] = -6
    omega_null_bc[i] = 0


for i in mesh.get_boundary_with_name('inflow'):
    psi_bc_value[i] = 3*y[i]**2 - 2*y[i]**3
    # vx[i] = v_in
    vx[i] = 6*y[i]-6*y[i]**2
    omega_bc_value[i] = 12*y[i]-6
    omega_null_bc[i] = 0

num_bc = len(mesh.boundary_nodes)


#Minv = sp.linalg.inv(M)
#Wz_old = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx))) * omega_null_bc
Wz_old = np.zeros(NN)


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


sp.random.seed(1)
break_flag = False
time_avg_loop = 0
iter = 0
for t in range(0, int(steps/vtk_steps)):
    for tt in range(vtk_steps):
        iter += 1
        time_start_loop = timer()
        print()
        print("Solving System " + str((float(iter)/(steps-1))*100) + "%")

        time_start_smooth = timer()
        vx_smooth, vy_smooth = m.weighted_smoothMesh(mesh.neighbour_nodes,
                                                     mesh.boundary_nodes,
                                                     x, y, dt)
        time_end_smooth = timer()
        print("Smoothing time: ", time_end_smooth - time_start_smooth)

        vx_w = sp.random.rand(NN) * 0.2*np.cos(2*np.pi*0.2*iter*dt)
        vy_w = sp.random.rand(NN) * 0.2*np.cos(2*np.pi*0.2*iter*dt)

        time_start_ale = timer()
        vx_mesh = p_lagrange * vx + p_smooth * vx_smooth + p_w * vx_w
        vy_mesh = p_lagrange * vy + p_smooth * vy_smooth + p_w * vy_w

        for i in mesh.boundary_nodes:
            vx_mesh[i] = 0
            vy_mesh[i] = 0

        x = x + vx_mesh * dt
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
        omega_bc = np.dot(Minv, (np.dot(Gx, vy) - np.dot(Gy, vx))) * omega_null_bc


        # Solve vorticity Transport
        LHS = M / dt + K / Re
        LHS_omega, omega_bc_RHS = apply_bc_dirichlet(LHS, mesh,
                                                     mesh.boundary_nodes,
                                                     omega_bc_value+omega_bc)

        F_omega = np.dot(M / dt, Wz_dep) + omega_bc_RHS
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
        err = Wz_new - Wz_old
        err = np.sqrt(np.dot(err,err))
        if err < 1e-8:
            break_flag = True
            break

        # Saving last step solution
        Psi_old = np.copy(Psi_new)
        Wz_old = np.copy(Wz_new)

        # Calculate Vx e Vy
    #    vx = sp.dot(Minv, sp.dot(Gy, Psi_new))
    #    vy = -1.0 * sp.dot(Minv, sp.dot(Gx, Psi_new))
        vx = sp.linalg.solve(M, np.dot(Gy, Psi_new))
        vy = -1.0 * sp.linalg.solve(M, np.dot(Gx, Psi_new))
        # Setting velocity BC
        for i in psi_dirichlet_nodes:
            vx[i] = 0
            vy[i] = 0
            if i in mesh.get_boundary_with_name('inflow'):
                vx[i] = 6*y[i]-6*y[i]**2


        time_end_loop = timer()
        time_avg_loop += time_end_loop - time_start_loop
        print("Iteration time: ", time_end_loop - time_start_loop)
        print("Estimated time to finish simulation: ",
                (time_end_loop - time_start_loop)* (steps-t+1) )


    # Save VTK
    vtk = io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, w_a,
                    vx_a, None, vx, vy, vxAle, vyAle)
    vtk.saveVTK(vtk_path,  's'+str(iter))
    vtk.saveVTK(vtk_path, 's-last')
    if break_flag: break

#-----------------  Loop End  --------------------
#-------------------------------------------------
time_end = timer()
print()
print("Iteration average time: ", time_avg_loop/(t+1))
print("Total runtime: ", time_end - time_start)
