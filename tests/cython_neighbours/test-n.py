import MESH as mm
import neighbours as nh
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import numpy as np


mesh = mm.Mesh('viv.msh')

it = 100

def run(n):
    t = np.zeros((n,3))
    for i in range(n):
        t1 = timer()
        a1, b1 = nh.py_neighbours(mesh.num_nodes, mesh.ien)
        t2 = timer()
        a2, b2 = nh.py_neighbours2(mesh.num_nodes, mesh.ien)
        t3 = timer()
        a3, b3 = nh.py_neighbours3(mesh.num_nodes, mesh.ien)
        t4 = timer()
        t[i][0] = t2-t1
        t[i][1] = t3-t2
        t[i][2] = t4-t3
    return t

t = run(it)
x = np.linspace(0, it+1, it)
plt.figure()
plt.plot(x, t[:,0], 'bo')
plt.plot(x, t[:,1], 'r--')
plt.plot(x, t[:,2], 'kx')
plt.show()
