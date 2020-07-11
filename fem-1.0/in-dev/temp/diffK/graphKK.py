# -*- coding: utf-8 -*-
import scipy as sp
import matplotlib.pyplot as plt
import Read as rd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
import os


data = rd.Read("KK-data")


#5kk
er5 = 0.0154449622735 / float(254)
l5 = 0.0953526623826 
#4kk
er4 = 0.0218600171775 / float(479)
l4 = 0.0694355572541
 
#3kk
er3 = 0.0295944409305 / float(851)
l3 = 0.0520936636879 
#2kk
er2 = 0.0387746828928 / float(1415)
l2 = 0.0403990787842 
#1kk
er1 = 0.0545908421987 / float(2790)
l1 = 0.0287705113166
 

error = [er5, er4, er3, er2, er1]
l = [l5, l4, l3, l2, l1]

linx = sp.linspace(0.0285,0.096,100)
liny = sp.zeros(100)
linx2 = sp.linspace(0.0285,0.096,100)
liny2 = sp.zeros(100)
for i in range(100):
	liny[i] = 0.03*linx[i]**2 
	liny2[i] = 0.0005*linx2[i] 

x = sp.linspace(0, 1, 101)

v_a = sp.zeros(1000)
x_a = sp.linspace(0, 1, 1000)
TA = 1/11.0
for i in range(1000):
	if x_a[i] > 0.5:
		v_a[i] = 2*(1-TA)*x_a[i] + 2*TA-1
	if x_a[i] <= 0.5:
		v_a[i] = 2*x_a[i]*TA 

v_a2 = sp.zeros(101)
for i in range(101):
	if x[i] > 0.5:
		v_a2[i] = 2*(1-TA)*x[i] + 2*TA-1
	if x[i] <= 0.5:
		v_a2[i] = 2*x[i]*TA 

erro_T = v_a2-data.T

print len(data.T), len(x), data.V




plt.figure(2)
plt.loglog(l, error, 'k-', marker='o')
plt.loglog(linx, liny, label =u'Quadrática')
plt.loglog(linx2, liny2, label ='Linear')
plt.grid('on')
plt.xlim(0.025, 0.1)
#plt.ylim(0.00007, 0.0012)
plt.title(u'Convergência')
plt.ylabel('Erro')
plt.xlabel('Tamanho do elemento')
plt.legend(loc='0')
plt.show()

exit()

plt.figure(3)
plt.title('Erro Absoluto')
plt.xlabel('y')
plt.ylabel('Erro')
plt.grid('on')
plt.xlim(-0.01, 1.01)
#plt.ylim(-0.01, 1.01)
plt.plot(x, erro_T,  'k.')


plt.figure(1)
plt.title('Temperatura')
plt.xlabel('y')
plt.ylabel('T')
plt.grid('on')
plt.xlim(-0.01, 1.01)
plt.ylim(-0.01, 1.01)
plt.plot(x, data.T,  'k.', label=u'Numérica')
plt.plot(x_a, v_a,  'b-', label=u'Analítica')
plt.legend(loc='0')




plt.show()


























