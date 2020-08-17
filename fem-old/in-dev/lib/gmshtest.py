import scipy as sp
import libmalha as lm
import InOut as io
import os
import random

cwd = os.getcwd()

def search2D(_neighbours, _ien, _x, _y, _xd, _yd, order='1'):
    px = _xd
    py = _yd
    end = 0
    j = random.randint(0,len(x)-1)
    dist = []
    result = -1
    while end == 0:
        for e in _neighbours[j]:

            x1 = _x[_ien[e, 0]]
            x2 = _x[_ien[e, 1]]
            x3 = _x[_ien[e, 2]]

            y1 = _y[_ien[e, 0]]
            y2 = _y[_ien[e, 1]]
            y3 = _y[_ien[e, 2]]

            A = sp.array([[x1, x2, x3],
                          [y1, y2, y3],
                          [1.0, 1.0, 1.0]])

            b = sp.array([px, py, 1.0])

            alpha = sp.linalg.solve(A, b)
            if sp.all(alpha >= 0.0) and sp.all(alpha <= 1.0):
                result = e
                end = 1
                break
            else:
                dist.append([_ien[e, 0], sp.sqrt((x1 - px) ** 2 + (y1 - py) ** 2)])
                dist.append([_ien[e, 1], sp.sqrt((x2 - px) ** 2 + (y2 - py) ** 2)])
                dist.append([_ien[e, 2], sp.sqrt((x3 - px) ** 2 + (y3 - py) ** 2)])
                end = 0

        if end == 1:
            break
        else:
            min_dist = min(dist, key=lambda k: k[1])
            last_node = j
            j = min_dist[0]

            if last_node == j and end == 0:
                print "point outside closer to ", last_node
                end = 1
                break
    print result
    return

x, y, ien = lm.lin2d(3,3,4,4)
ien = lm.ret_tri(ien)
v = sp.zeros((len(x),2))


nei = lm.neighbourElements(len(x), ien)

vtk = io.InOut(x, y, ien, len(x), len(ien), None,   v[:, 0], v[:, 1])
vtk.saveVTK(cwd, "test")

px = random.randint(-10,60)*0.1
py = random.randint(-10,50)*0.1

search2D(nei, ien, x, y, px, py)

print px, py
