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


erro_eu = sp.array([1.90884940e-03, 1.76660643e-04, 4.67488783e-05, 1.74977724e-05, 6.23026825e-06, 2.44087891e-06])
logerro_eu = sp.zeros(6)
l = sp.array([0.8, 0.3, 0.15, 0.1, 0.08, 0.05])

print  erro_eu

linx = sp.linspace(0.05,0.8,100)
liny = sp.zeros(100)
linx2 = sp.linspace(0.05,0.8,100)
liny2 = sp.zeros(100)
for i in range(100):
	liny[i] = 0.01*linx[i]**2 
	liny2[i] = 0.00001*linx2[i] 


plt.figure(3)
plt.loglog(l, erro_eu, 'k-', marker='o')
plt.loglog(linx, liny, label ='Quadratic')
plt.loglog(linx2, liny2, label ='Linear')
plt.grid('on')
plt.title('Convergence')
plt.ylabel('Error')
plt.xlabel(r'$l$')
plt.legend(loc='0')
plt.show()


