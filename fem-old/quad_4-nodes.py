## =================================================================== ##
#  this is file quad_4-nodes.py, created at 30-Apr-2020                #
#  maintained by Gustavo Rabello dos Anjos                              #
#  e-mail: gustavo.rabello@gmail.com                                    #
## =================================================================== ##

from sympy import *
init_printing(use_unicode=False, wrap_line=False)

x = Symbol('x')
y = Symbol('y')
a = Symbol('a')
b = Symbol('b')
kx = Symbol('kx')
ky = Symbol('ky')

phi = [0,0,0,0]
phi[0] = (1-x/(2*b)) * (1-y/(2*a))
phi[1] = x/(2*b) * (1-y/(2*a))
phi[2] = x*y/(4*a*b)
phi[3] = y/(2*a) * (1-x/(2*b))

dphidx = [0,0,0,0]
dphidx[0] = diff(phi[0],x)
dphidx[1] = diff(phi[1],x)
dphidx[2] = diff(phi[2],x)
dphidx[3] = diff(phi[3],x)

dphidy = [0,0,0,0]
dphidy[0] = diff(phi[0],y)
dphidy[1] = diff(phi[1],y)
dphidy[2] = diff(phi[2],y)
dphidy[3] = diff(phi[3],y)

# Mass matrix
m = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
for i in range(0,4):
 for j in range(0,4):
  g=integrate(phi[i]*phi[j],(x,0,2*b))
  f=integrate(g,(y,0,2*a))
  m[i,j] = f

# Gradient_x matrix
gx = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
for i in range(0,4):
 for j in range(0,4):
  g=integrate(phi[i]*dphidx[j],(x,0,2*b))
  f=integrate(g,(y,0,2*a))
  gx[i,j] = f

# Gradient_y matrix
gy = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
for i in range(0,4):
 for j in range(0,4):
  g=integrate(phi[i]*dphidy[j],(x,0,2*b))
  f=integrate(g,(y,0,2*a))
  gy[i,j] = f

# Stiffness matrix
kxx = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
for i in range(0,4):
 for j in range(0,4):
  g=integrate(dphidx[i]*dphidx[j],(x,0,2*b))
  f=integrate(g,(y,0,2*a))
  kxx[i,j] = f

kyy = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
for i in range(0,4):
 for j in range(0,4):
  g=integrate(dphidy[i]*dphidy[j],(x,0,2*b))
  f=integrate(g,(y,0,2*a))
  kyy[i,j] = f

kxy = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
for i in range(0,4):
 for j in range(0,4):
  g=integrate(dphidx[i]*dphidy[j],(x,0,2*b))
  f=integrate(g,(y,0,2*a))
  kxy[i,j] = f


# mass matrix
print ("mass:")
print ("2ab/18 *")
pprint (18/(2*a*b)*m)
print ("")

# stiffness kxx matrix
print ("a/6b *")
print ("stiffness kxx:")
pprint ((6*b)/(a) * kxx)
print ("")

# stiffness kyy matrix
print ("b/6a *")
print ("stiffness kyy:")
pprint ((6*a)/(b) * kyy)
print ("")

# stiffness kxy matrix
print ("stiffness kxy:")
print ("1/4 *")
pprint (4 * kxy)
print ("")

# stiffness gx matrix
print ("a/6 *")
print ("gradient Gx:")
pprint ((6)/(a) * gx)
print ("")

# stiffness gy matrix
print ("gradient Gy:")
print ("b/6 *")
pprint ((6)/(b) * gy)
print ("")


