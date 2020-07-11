# -*- coding: utf-8 -*-
import scipy as sp
from scipy import special
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib import cm


# -------------------------------------------------------
#  Mesh
# -------------------------------------------------------

L = 5		# cylinder length
R = 1		# cylinder radius
num_x = 10 	# number of divisions in z-direction
num_y = 10	# number of divisions in r-direction

dx = float(L) / num_x
dy = float(R) / num_y
x = sp.zeros((num_x + 1) * (num_y + 1))		# z-direction point coordinates
y = sp.zeros((num_x + 1) * (num_y + 1))		# r-direction point coordinates
k = 0
for j in range(0, (num_x + 1) * (num_y + 1), num_x + 1):
    for i in range(0, num_x + 1):
        x[i + j] = i * dx
for i in range(0, (num_x + 1) * (num_y + 1), num_x + 1):
    for j in range(0, num_x + 1):
        y[i + j] = k * dy
    k += 1

nodes = len(x)

# -------------------------------------------------------
#  Main
# -------------------------------------------------------

q = 1		# Heat Generation
k = 1		# Thermal Conductivity

maxsum = 100		# Number of infinite sum terms

bc = 1		# Constant boundary value

# Exemple Boundary function:
def boundary(_z):	
    return (1./float(L)) * _z


def gz_sin(_z, _ln, _bcFunction = bc):
    return (_bcFunction - (q*L**2)*(1./(2.*k)) * ((_z/L) - (_z/L)**2)) * sp.sin(_ln * _z)

temperature = sp.zeros(nodes)
for i in range(nodes):
    print i, " / ", nodes
    z = x[i]
    r = y[i]
    for n in range(1, maxsum):
        ln = n*(sp.pi/5.)
        i0 = sp.special.iv(0, ln)
        ir = sp.special.iv(0, r*ln)
        c1 = sp.integrate.quadrature(lambda e: gz_sin(e, ln), 0, L, maxiter=100)[0]      # constant boundary bc
        # c1 = sp.integrate.quadrature(lambda e: gz_sin(e, ln, boundary(e)), 0, L, maxiter=100)[0]  # Boundary function
        c2 = 0.5 * L - sp.sin(2*ln*L) / (4*ln)
        cn = c1 / (i0 * c2)
        temperature[i] += cn * ir * sp.sin(z*ln)
    temperature[i] += (q*L**2)*(1./(2.*k)) * ((z/L) - (z/L)**2)

plt.scatter(x, y, c=temperature, cmap=cm.jet, marker="s", linewidth=10)
plt.colorbar()
plt.show()

