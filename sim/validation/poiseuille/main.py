import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# A = area_t / NE
# r = np.sqrt( A / (3*np.sqrt(3)))
# l = np.sqrt( 4 * (5/N) / (np.sqrt(3)))

# h = [0.11724519171678374, 0.07352483523889468, 0.05045429721577769,
     # 0.02595852176059649, 0.01572380485828159]

h = [0.40614925799324625, 0.25469750050379236, 0.1747788124758158,
     0.08992295715747085, 0.054468857805684126]

# h = [70, 178, 378, 1428, 3892]

# ---------------------------------------
# Analytical Solution
# ---------------------------------------

def get_analytic(_x, _y):
    _num = len(_x)
    Wz_a = np.zeros(_num, dtype="float64")
    Psi_a = np.zeros(_num, dtype="float64")
    vx_a = np.zeros(_num, dtype="float64")

    for i in range(_num):
        vx_a[i] = 6 * _y[i] * (1 - _y[i])
        Psi_a[i] = 3*_y[i]**2 - 2*_y[i]**3
        Wz_a[i] = 12*_y[i] - 6

    return Wz_a, Psi_a, vx_a


error = np.zeros((3,5), dtype="float64")

for t in range(5):
    s = pd.read_csv('{}.csv'.format(t+1))
    x = list(s['Points:0'])
    y = list(s['Points:1'])
    n = len(x)
    wa, pa, va = get_analytic(x,y)
    error[0,t] = np.sqrt(np.square(list(s['Omega'])- wa).sum()) / n
    error[1,t] = np.sqrt(np.square(list(s['Psi'])- pa).sum()) /n
    error[2,t] = np.sqrt(np.square(list(s['Velocity:0'])- va).sum()) /n


linx = np.linspace(0.05,0.8,100)
liny = np.zeros(100)
linx2 = np.linspace(0.05,0.8,100)
liny2 = np.zeros(100)
for i in range(100):
	liny[i] = 0.8*linx[i]**2
	liny2[i] = 2*linx2[i]

plt.figure(1)
plt.loglog(h, error[2], 'k-', marker='o', label = 'Numeric')
plt.loglog(linx, liny, 'k--', label ='Quadratic')
plt.loglog(linx2, liny2, 'k-',label ='Linear')
plt.grid('on')
plt.title('Vorticity Convergence Rate')
plt.ylabel('Error')
plt.xlabel(r'$h$')
plt.legend(loc='0')
plt.show()
