# -*- coding: utf-8 -*-
import numpy as np
from libmalha import lin2d, ret_tri
from timeit import default_timer as timer
import matplotlib.pyplot as plt
from assembly import fem_matrix as fem_py
from cyassembly import fem_matrix as fem_cy


L = 10.0


def test(N):

    L = 10.0
    x,y,ien = lin2d(N,L,N,L)
    ien = ret_tri(ien)

    nodes = len(x)
    elem = len(ien)

    # fem_matrix(_x, _y, _numele, _numnode, _ien)

    t1 = timer()
    k,m,gx,gy = fem_py(x, y, elem, nodes, ien)
    t2 = timer()
    k,m,gx,gy = fem_cy(x, y, elem, nodes, ien)
    t3 = timer()
    return t2-t1, t3-t2, elem

size = 70
pyt = np.zeros(size)
cyt = np.zeros(size)
ns = np.zeros(size)
for n in range(1,size+1):
    print(n)
    pyt[n-1], cyt[n-1], ele = test(n)
    ns[n-1] = ele

plt.figure(1)
plt.title('FEM Matrix Assembly timing')
plt.semilogx(ns, pyt, 'k-', label = 'Python')
plt.semilogx(ns, cyt, 'b-', label = 'Cython')
plt.ylabel('Time [s]')
plt.xlabel('Number of elements')
plt.legend(loc=0)
plt.show()
