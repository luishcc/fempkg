# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
import os

cwd = os.getcwd()


# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

arquivo = "1kk"

malha_total = Gm.GMesh(arquivo+".msh")

x_total = malha_total.X
y_total = malha_total.Y
ien_total = malha_total.IEN
nodes_total = len(x_total)
num_ele_total = len(ien_total)




# -------------------------------------------------------
#     Defining Parameters
# -------------------------------------------------------

dt = 0.005
tempo = 400

alfa_solid = 1.0

cp_fld = 1.0
rho_fld = 1.0
kt_fld = 0.1

alfa_fld = kt_fld / (rho_fld*cp_fld)


# ---------------------------------------
# Wz, Psi, T e velocidade inicial
# ---------------------------------------


Temp_old = sp.ones(nodes_total, dtype="float64")
for i in range(len(malha_total.dirichlet_points)):
    index = int(malha_total.dirichlet_points[i][0]-1)
    value = malha_total.dirichlet_points[i][1]
    Temp_old[index] = value


# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_local = sp.zeros((3, 3), dtype="float64")
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    #gx_local = sp.zeros((3, 3), dtype="float64")
    #gy_local = sp.zeros((3, 3), dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    #gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    #gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):
        alfa = 0
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]
	    
	if (sp.amax(yy) > 0.5) and (sp.amin(yy) >= 0.5): 
            alfa = alfa_fld
	else:
	    alfa = alfa_solid

	print sp.amax(yy), alfa

       # alfa = alfa * (1/3.0)
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
                #gx_local[i,j] = b[j] * (1/6.)
                #gy_local[i,j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]*alfa
                m_global[i_global, j_global] += m_local[i_local, j_local]* (Area/12.)
                #gx_global[i_global, j_global] += gx_local[i_local, j_local]
                #gy_global[i_global, j_global] += gy_local[i_local, j_local]


    return  k_global, m_global#, gx_global, gy_global

# --------------------------------------
# Total Matrices (Heat equation)


Alfa_total = sp.zeros(nodes_total)
for i in range(nodes_total):
    if y_total[i] > 0.5:
        Alfa_total[i] = alfa_fld
    if y_total[i] < 0.5:
        Alfa_total[i] = alfa_solid
    if y_total[i] == 0.5:
        Alfa_total[i] = (2*alfa_solid*alfa_fld) / (alfa_fld+alfa_solid)

K, M = fem_matrix(x_total, y_total, num_ele_total, nodes_total, ien_total)

Mdt = M/dt

# --------------------------------------


# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------

dirichlet_len_total = len(malha_total.dirichlet_points)


T_a = sp.zeros(nodes_total)
Tan = 1/11.0
for i in range(nodes_total):
    if y_total[i] > 0.5:
        T_a[i] = 2 * (1 - Tan) * y_total[i] + 2 * Tan - 1
    else:
        T_a[i] = 2 * y_total[i] * Tan


# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

for t in range(0, tempo):
    print "Solving System " + str((float(t)/(tempo))*100) + "%"


    LHS_T = Mdt + K
    LHS_temp = sp.copy(LHS_T)
    cctemp = sp.zeros(nodes_total)
    for i in range(dirichlet_len_total):
        index = int(malha_total.dirichlet_points[i][0]) - 1
        value = malha_total.dirichlet_points[i][1]
        for j in range(nodes_total):
            cctemp[j] -= value * LHS_T[j, index]
            if j != index:
                LHS_temp[index, j] = 0
                LHS_temp[j, index] = 0
            else:
                LHS_temp[index, j] = 1

    F_temp = sp.dot(Mdt, Temp_old) + cctemp
    for i in range(dirichlet_len_total):
        index = int(malha_total.dirichlet_points[i][0]) - 1
        F_temp[index] = malha_total.dirichlet_points[i][1]

    Temp_new = sp.linalg.solve(LHS_temp, F_temp)

    # Salvar VTK
    vtk_t = Io.InOut(x_total, y_total, ien_total, nodes_total, num_ele_total, Temp_old, T_a, None, None, None, None, None)
    vtk_t.saveVTK(cwd+"/results", arquivo + str(t+1))

    Temp_old = sp.copy(Temp_new)

# ----------------- Fim de Loop -------------------
# -------------------------------------------------


error = Temp_new - T_a
l2norm = 0.0
for i in range(nodes_total):
    l2norm += error[i]**2
l2norm = sp.sqrt(l2norm)

Ele_area = 1/float(nodes_total)
char_l = sp.sqrt(Ele_area * (4/3.) * sp.sqrt(3))

print l2norm, char_l



fig = plt.figure(1)
ax = fig.gca(projection='3d')
plt.title(u'Erro Absoluto')
plt.axis('on')
#ax.set_xlim3d(0, 3)
#ax.set_ylim3d(0, 2)
#ax.set_zlim3d(0, 1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('E')
ax.view_init(24, -135)
ax.plot_trisurf(x_total, y_total, sp.absolute(error), cmap='jet')
plt.savefig('t_kk_abs.png')



# ----------------------------------
