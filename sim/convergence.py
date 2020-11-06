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


psi_top=1

cwd = os.getcwd()

num_case=5

msh_file = ["poiseuille-{}".format(str(i+1)) for i in range(num_case)]
sim_case = 'convergence'

results_path = rt_save.make_dir(sim_case)

print('-------------------------------------------------------------')
print()
print('  Starting Simulation: {}'.format(sim_case))
print('  Saving Results in: {}'.format(results_path))
print()
print('-------------------------------------------------------------')

def solve(mesh):

    x = mesh.x
    y = mesh.y
    ien = mesh.ien
    NN = mesh.num_nodes
    NE = mesh.num_elem

    psi = np.zeros(NN, dtype="float64")
    vx = np.zeros(NN, dtype="float64")

    vx_a = np.zeros(NN, dtype="float64")
    w_a = np.zeros(NN, dtype="float64")
    psi_a = np.zeros(NN, dtype="float64")
    for i in range(NN):
        vx_a[i] = 6*y[i]*(1-y[i])
        w_a[i] = 12*y[i] - 6
        psi_a[i] = 3*y[i]**2 - 2*y[i]**3

    K, M, Gx, Gy = fem_matrix(mesh)

    psi_dirichlet_nodes = mesh.set_dirichlet_nodes(['top', 'bot',
                                                    'inflow'])
    psi_bc_value = np.zeros(NN)

    for i in mesh.get_boundary_with_name('top'):
        psi_bc_value[i] = 1
        vx[i] = 0

    for i in mesh.get_boundary_with_name('bot'):
        psi_bc_value[i] = 0
        vx[i] = 0

    for i in mesh.get_boundary_with_name('inflow'):
        psi_bc_value[i] = psi_a[i]
        vx[i] = vx_a[i]

    num_bc = len(mesh.boundary_nodes)

    K_psi, psi_bc_RHS = apply_bc_dirichlet(K, mesh, psi_dirichlet_nodes,
                                            psi_bc_value)
    F_psi = np.dot(M, w_a) + psi_bc_RHS
    for i in psi_dirichlet_nodes:
        F_psi[i] = psi_bc_value[i]
    psi = sp.linalg.solve(K_psi, F_psi)

    # LHS_vx, BC_vx = apply_bc_dirichlet(M, mesh, psi_dirichlet_nodes, vx_a)
    # F_vx = np.dot(Gy, psi) + BC_vx
    # for i in psi_dirichlet_nodes:
    #     F_vx[i] = vx_a[i]
    # vx = sp.linalg.solve(LHS_vx, F_vx)
    vx = sp.linalg.solve(M, np.dot(Gy, psi))

    wz = sp.linalg.solve(M, (- np.dot(Gy, vx)))

    err_vx = vx - vx_a
    err_vx = np.sqrt(np.dot(err_vx,err_vx))

    err_psi = psi - psi_a
    err_psi = np.sqrt(np.dot(err_psi, err_psi))

    err_w = wz - w_a
    err_w = np.sqrt(np.dot(err_w, err_w))

    print(err_vx, err_psi, np.sqrt( (4 * (5./NE)) / np.sqrt(3)))
    return err_vx/NE, err_psi/NE, np.sqrt( (4 * (5./NE)) / np.sqrt(3) ), err_w/NE
    # return err_vx/NE, err_psi/NE, np.sqrt(5./NE), err_w/NE
    # return err_vx/NE, err_psi/NE, 5./NE, err_w/NE
    # return err_vx/NE, err_psi/NE, np.sqrt( (5./NE) /(3* np.sqrt(3)) ), err_w/NE




e_v = np.zeros(num_case)
e_p = np.zeros(num_case)
e_w = np.zeros(num_case)
ar = np.zeros(num_case)

for t in range(num_case):
    mm = m.Mesh("mesh/{}/{}.msh".format(sim_case,msh_file[t]))
    e_v[t], e_p[t], ar[t], e_w[t] = solve(mm)

linx = np.linspace(0.03,0.5,100)
liny = np.zeros(100)
linx2 = np.linspace(0.03,0.5,100)
liny2 = np.zeros(100)
for i in range(100):
	liny[i] = 0.3*linx[i]**2
	liny2[i] = 0.8*linx2[i]

# A = area_t / NE
# r = np.sqrt( A / (3*np.sqrt(3)))
# l = np.sqrt( 4 * (5/N) / (np.sqrt(3)))

# h = [0.11724519171678374, 0.07352483523889468, 0.05045429721577769,
     # 0.02595852176059649, 0.01572380485828159]

h = [0.40614925799324625, 0.25469750050379236, 0.1747788124758158,
     0.08992295715747085, 0.054468857805684126]

#ar = h[0:num_case]

# h = [70, 178, 378, 1428, 3892]


import matplotlib.pyplot as plt

plt.figure(1)
plt.loglog(ar, e_v, 'k-', marker='o', label = 'Numeric')
plt.loglog(linx, 0.5*liny*1, 'k--', label ='Quadratic')
plt.loglog(linx2, liny2*1, 'k-',label ='Linear')
plt.grid('on')
plt.title('Velocity Convergence Rate')
plt.ylabel('Error')
plt.xlabel(r'$h$')
plt.legend()

plt.figure(2)
plt.loglog(ar, e_p, 'k-', marker='o', label = 'Numeric')
plt.loglog(linx, liny/100, 'k--', label ='Quadratic')
plt.loglog(linx2, liny2/100, 'k-',label ='Linear')
plt.grid('on')
plt.title('Stream Function Convergence Rate')
plt.ylabel('Error')
plt.xlabel(r'$h$')
plt.legend()


plt.figure(3)
plt.loglog(ar, e_w, 'k-', marker='o', label = 'Numeric')
plt.loglog(linx, liny*10, 'k--', label ='Quadratic')
plt.loglog(linx2, liny2*10, 'k-',label ='Linear')
plt.grid('on')
plt.title('Vorticity Convergence Rate')
plt.ylabel('Error')
plt.xlabel(r'$h$')
plt.legend()

plt.show()
