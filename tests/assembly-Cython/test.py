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
plt.ylabel('Wall Time [s]')
plt.xlabel('Number of elements')
plt.legend(loc=0)
plt.show()

#-------------------------------------------------------------
#-------------------------------------------------------------
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 0.1 #set the value globally
mpl.rcParams['lines.dashed_pattern'] = [6, 6]
mpl.rcParams['lines.scale_dashes'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import warnings
warnings.simplefilter('ignore')

from matplotlib import pyplot as plt
from scipy.sparse.linalg.dsolve import linsolve



fig1, ax = plt.subplots(1,1,figsize=(3.5,3), frameon=True)
ax.locator_params(axis = 'y',nbins=10)
ax.locator_params(axis = 'x',nbins=10)
ax.tick_params(direction='in', width='0.1')
plt.tick_params(axis='both', labelsize=12, pad=1)
plt.xlabel(r'N', fontsize=12,labelpad=2)
plt.ylabel(r'Wall time $[s]$',style='italic', fontsize=12, rotation='vertical',
           labelpad=3)
ax.grid(color='grey', linestyle='--', linewidth=0.1,
        drawstyle='steps')
#plt.plot(A[:,1],A[:,0],'r-', linewidth=0.75,label='python')
plt.plot(ns, pyt,'b-', linewidth=0.75,label='Python')
plt.plot(ns, cyt,'ko', markersize=1.5,label='Cython')

plt.legend(title='',loc=1, borderaxespad=0., fontsize=12)
legend = plt.legend(loc=2, fontsize=10, facecolor='white',frameon=True,
                   fancybox=False, framealpha=1, shadow=False, borderpad=0.2,
                    labelspacing=0.5)
legend.get_frame().set_linewidth(0.1)
legend.get_frame().set_edgecolor("black")
plt.title('Assembly Comparison')
ax.set_xscale('log')
plt.show()
fig1.savefig('fem_cy_v_py.pdf',
             format='pdf',dpi=300, bbox_inches='tight',
             transparent=True)
