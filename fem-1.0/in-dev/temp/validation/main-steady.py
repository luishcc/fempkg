# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import GMesh as gm
import numpy as np
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

# Definição da malha
arquivo_pre = "validation"

erro_eu = sp.zeros(6)
logerro_eu = sp.zeros(6)
l = sp.array([0.8, 0.3, 0.15, 0.1, 0.08, 0.05])


def fem_matrix(_x, _y, _numele, _numnode, _ien, _alfa):
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
            alfa += _alfa[_ien[elem, i]]

        alfa = alfa * (1/3.0)
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


    return  k_global, m_global #, gx_global, gy_global


for arq in range(1,7):

	arquivo = arquivo_pre + str(arq)

	print "Reading "+arquivo+".msh file"
	malha = gm.GMesh(arquivo+".msh")

	x = malha.X
	y = malha.Y
	ien = malha.IEN

	# Parametros
	nodes = len(x)
	elenum = len(ien)
	dt = 0.1
	tempo = 100
	teta = 1        # [0, 1] --> explicito - implicito


	Alfa = sp.ones(nodes)

	# Montagem de matrizes

	print "Assembling Matrices"


	#K, M, Gx, Gy = fem_matrix(x, y, elenum, nodes, ien)

	K, M = fem_matrix(x, y, elenum, nodes, ien, Alfa)


	    # Implementando condições de contorno de DIRICHLET e INICIAL

	print "Defining Boundary and Initial conditions"


	quantcc = len(malha.dirichlet_points)
	cc = sp.zeros(nodes)


	T_old = sp.zeros(nodes)
	LHS = K
	KK = sp.copy(LHS)

	for i in range(quantcc):
	    index = int(malha.dirichlet_points[i][0]-1)
	    value = malha.dirichlet_points[i][1]

	    if x[index] == 3.0 and y[index] >= 1.0:
		value = 2-y[index]

	    if x[index] == 3.0 and y[index] < 1.0:
		value = y[index]

	    T_old[index] = malha.dirichlet_points[i, 1]
	    for j in range(nodes):
		cc[j] -= value * LHS[j, index]
		if j != index:
		    KK[index, j] = 0
		    KK[j, index] = 0
		else:
		    KK[index, j] = 1


	    # Solucao do sistema

	T_new = sp.zeros(nodes)


	T_a = sp.zeros(nodes, dtype="float128")
	pi = np.pi
	for i in range(nodes):
	    T_a[i] = 0
	    for m in range(1,100):
		#print m
		T_a[i] += (8.0/pi**2) * np.sin(m*0.5*pi) * (1./(m**2 * np.sinh(1.5*m*pi))) * np.sinh(m*0.5*pi*x[i]) * np.sin(m*(1./2.)*pi*y[i])


	B = cc
	for j in range(quantcc):
		index = int(malha.dirichlet_points[j][0]) - 1
		value = malha.dirichlet_points[j][1]
		if x[index] == 3.0 and y[index] >= 1.0:
		    value = 2 - y[index]
		if x[index] == 3.0 and y[index] < 1.0:
		    value = y[index]
		B[index] = value
	T_new = linalg.solve(KK, B)


	T_old = sp.copy(T_new)



	erro_abs = sp.zeros(nodes)
	erro_eu[arq-1] = 0.
	for i in range(nodes):
		erro_abs[i] = (T_old[i]-T_a[i])
		erro_eu[arq-1] += erro_abs[i]**2

	erro_eu[arq-1] = sp.sqrt(erro_eu[arq-1])/nodes
	logerro_eu[arq-1] = sp.log10(erro_eu[arq-1])

print  erro_eu

vtk = io.InOut(x, y, ien, len(x), len(ien), T_old, T_a, None)
vtk.saveVTK(cwd , arquivo + str(i))

linx = sp.linspace(0.05,0.8,100)
liny = sp.zeros(100)
linx2 = sp.linspace(0.05,0.8,100)
liny2 = sp.zeros(100)
for i in range(100):
	liny[i] = 0.01*linx[i]**2 
	liny2[i] = 0.00001*linx2[i] 





fig = plt.figure(4)
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
ax.plot_trisurf(x, y, erro_abs, cmap='jet')
plt.savefig(arquivo+'_abs.png')
plt.show()



plt.figure(3)
plt.loglog(l, erro_eu, 'k-', marker='o')
plt.loglog(linx, liny, label =u'Quadrática')
plt.loglog(linx2, liny2, label ='Linear')
plt.grid('on')
plt.title(u'Convergência')
plt.ylabel('Erro')
plt.xlabel('Tamanho do elemento$')
plt.legend(loc='0')
plt.savefig(arquivo+'_conv.png')


exit()
fig = plt.figure(1)
ax = fig.gca(projection='3d')
plt.title(u'Solução Numérica')
plt.axis('on')
ax.set_xlim3d(0, 3)
ax.set_ylim3d(0, 2)
ax.set_zlim3d(0, 1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u')
ax.view_init(24, -135)
ax.plot_trisurf(x, y, T_old, cmap='jet')
plt.savefig(arquivo+'_num.png')


fig = plt.figure(2)
ax2 = fig.gca(projection='3d')
plt.axis('on')
ax2.set_xlim3d(0, 3)
ax2.set_ylim3d(0, 2)
ax2.set_zlim3d(0, 1)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('u')
plt.title(u'Solução Analítica')
ax2.view_init(24, -135)
ax2.plot_trisurf(x, y, T_a, cmap='jet')
plt.savefig(arquivo+'_alyt.png')





