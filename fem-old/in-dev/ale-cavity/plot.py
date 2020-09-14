import numpy as np
import matplotlib.pyplot as plt

# Ghia

y_a = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
x_a = [1.0000, 1.0000, 0.9688, 0.9609, 0.9531, 0.9060, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000 ]
vx_a = [ 1.00000,  0.84123,  0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717,  0.00000 ]
vy_a = [ 0.00000, 0.00000, -0.0739, -0.0886, -0.1031, -0.1691, -0.2244, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000 ]


# Marchi

y_b = [ 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000, 0.4375, 0.3750, 0.3125, 0.2500, 0.1875, 0.1250, 0.0625, 0.0000 ]
x_b = [ 1.0000, 0.9375, 0.8750, 0.8125, 0.7500, 0.6875, 0.6250, 0.5625, 0.5000, 0.4375, 0.3750, 0.3125, 0.2500, 0.1875, 0.1250, 0.0625, 0.0000 ]
vx_b = [ 1.00000, 0.59746, 0.31055, 0.14042, 0.02787, -0.0602, -0.1312, -0.18208, -0.20914, -0.21296, -0.19847, -0.17271, -0.14193, -0.10981, -0.07712, -0.04197, 0.00000 ]
vy_b =[ 0.00000, -0.1233,  -0.2186,  -0.2537,  -0.2278, -0.1630, -0.0840, -0.00774, 0.05753, 0.10877, 0.14573, 0.16913, 0.17924, 0.17434, 0.14924,  0.09480, 0.00000 ]


# Sim

import pandas as pd

horizontal = pd.read_csv('horizontal.csv')
vertical = pd.read_csv('vertical.csv')


vx_s = list(vertical['Velocity:0'])
vy_s = list(horizontal['Velocity:1'])
p_s = np.linspace(0,1,len(vx_s))



plt.figure(1)
plt.plot(x_a, vy_a, 'ko', label='Ghia et al.')
plt.plot(x_b, vy_b, 'kx', label='Marchi et al.')
plt.plot(p_s, vy_s, label='ALE' )
plt.title('Plot Over Horizontal Line')
plt.xlabel('x')
plt.ylabel('Vertical Velocity')
plt.legend()

plt.figure(2)
plt.plot(y_a, vx_a, 'ko', label='Ghia et al.')
plt.plot(y_b, vx_b, 'kx', label='Marchi et al.')
plt.plot(p_s, vx_s, label='ALE')
plt.title('Plot Over Vertical Line')
plt.xlabel('y')
plt.ylabel('Horizontal Velocity')
plt.legend()
plt.show()



