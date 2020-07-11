# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
import os
import semiLagrangean as sl
import matplotlib.pyplot as plt
cwd = os.getcwd()

mesh_file = "test"
# mesh_file = "test-s"
savename = mesh_file


mesh = gm.GMesh("mesh/" + mesh_file + ".msh")

num_nodes = len(mesh.X)
num_ele = len(mesh.IEN)

neighbour_ele, neighbour_nodes = sl.neighbourElements2(num_nodes, mesh.IEN)

# print mesh.Boundary_Nodes

def smoothMesh(_neighbour_nodes, _mesh, _x, _y, _dt, delT = 1):

    xx = sp.copy(_x)
    yy = sp.copy(_y)

    for i in range(len(_neighbour_nodes)):

        flag = 0
        for k in _mesh.Boundary_Nodes:
            if i == k:
                flag = 1
                break

        if flag == 1:
            continue

        vertex_position = sp.array([xx[i], yy[i]])
        nghN = _neighbour_nodes[i]
        num_nghb = len(nghN)
        distance_vectors = sp.zeros((num_nghb, 2))
        displacement_vector = sp.zeros(2)
        for j in range(num_nghb):
            index = nghN[j]
            distance_vectors[j][0] = xx[index] - vertex_position[0]
            distance_vectors[j][1] = yy[index] - vertex_position[1]
            displacement_vector += distance_vectors[j]
        displacement_vector *= (1./num_nghb)
        displacement_velocity = displacement_vector / delT


        #new_position = vertex_position + displacement_vector
        new_position = vertex_position + displacement_velocity * _dt

        xx[i] = new_position[0]
        yy[i] = new_position[1]


      #  print nghN, "\n", distance_vectors, "\n", displacement_vector, "\n", new_position

    return xx, yy

def weightedSmoothMesh(_neighbour_nodes, _mesh, _x, _y, _dt, delT = 1):

    xx = sp.copy(_x)
    yy = sp.copy(_y)

    for i in range(len(_neighbour_nodes)):

        flag = 0
        for k in _mesh.Boundary_Nodes:
            if i == k:
                flag = 1
                break
        if flag == 1:
            continue

        vertex_position = sp.array([xx[i], yy[i]])
        nghN = _neighbour_nodes[i]
        num_nghb = len(nghN)
        neighbour_positions = sp.zeros((num_nghb, 2))
        displacement_vector = sp.zeros(2)
        sumWeight = 0
        for j in range(num_nghb):
            index = nghN[j]
            neighbour_positions[j][0] = xx[index]
            neighbour_positions[j][1] = yy[index]
            distance = sp.sqrt((xx[index] - vertex_position[0])**2 +
                               (yy[index] - vertex_position[1])**2)
#            distance = 1

            # displacement_vector += (1./distance) * neighbour_positions[j]
            # sumWeight += (1./distance)

            displacement_vector += (distance) * neighbour_positions[j]
            sumWeight += (distance)

        displacement_vector = (displacement_vector/sumWeight) - vertex_position

        displacement_velocity = displacement_vector / delT

        new_position = vertex_position + displacement_velocity * _dt

        xx[i] = new_position[0]
        yy[i] = new_position[1]

      #  print nghN, "\n", distance_vectors, "\n", displacement_vector, "\n", new_position

    return xx, yy

def distance(_a,_b):
    size = len(_a)
    sum = 0
    for i in range(size):
        sum += (_b[i] - _a[i])**2
    return sp.sqrt(sum)

def triArea(_a, _b, _c):
    sa = distance(_b, _c)
    sb = distance(_a, _c)
    sc = distance(_b, _a)
    s = 0.5 * (sa + sb + sc)
    return sp.sqrt(s * (s-sa) * (s-sb) * (s-sc))

def triAreaRelation(_a, _b, _c):
    sa = distance(_b, _c)
    sb = distance(_a, _c)
    sc = distance(_b, _a)
    s = 0.5 * (sa + sb + sc)
    area = sp.sqrt(s * (s-sa) * (s-sb) * (s-sc))
    return (3 * area * sp.sqrt(3)) / s**2

def checkMesh(_ien, _x, _y):
    sum = 0
    for e in range(len(_ien)):
        v1 = [_x[_ien[e][0]], _y[_ien[e][0]]]
        v2 = [_x[_ien[e][1]], _y[_ien[e][1]]]
        v3 = [_x[_ien[e][2]], _y[_ien[e][2]]]
        sum += triAreaRelation(v1, v2, v3)
    return sum / len(_ien)


x2 = sp.copy(mesh.X)
y2 = sp.copy(mesh.Y)

x3 = sp.copy(mesh.X)
y3 = sp.copy(mesh.Y)

vtk_t = IO.InOut(x2, y2, mesh.IEN, num_nodes, num_ele, None, None, None, None, None, None, None)
vtk_t.saveVTK(cwd + "/results", savename + str(0))


iter = 2000
t2 = sp.zeros(iter)
t3 = sp.zeros(iter)
for t in range(1, iter):
    print "Iteration ", t, "/", iter
    x2, y2 = smoothMesh(neighbour_nodes, mesh, x2, y2, 0.001, delT=0.1)
    x3, y3 = weightedSmoothMesh(neighbour_nodes, mesh, x3, y3, 0.001, delT=0.1)

    t2[t] = checkMesh(mesh.IEN, x2, y2)
    t3[t] = checkMesh(mesh.IEN, x3, y3)

    vtk_t = IO.InOut(x2, y2, mesh.IEN, num_nodes, num_ele, None, None, None, None, None, None, None)
    vtk_t.saveVTK(cwd + "/results", savename + str(t))

t1 = checkMesh(mesh.IEN, mesh.X, mesh.Y)
t2[0] = t1
t3[0] = t1
print t1

plt.figure(1)
plt.title('Mesh Smoothing Technique')
plt.plot(sp.linspace(0, iter, iter), t2, 'k', label='Laplacian')
plt.plot(sp.linspace(0, iter, iter), t3, 'b', label='Weighted Laplacian')
plt.legend(loc='0')
plt.xlabel('iterations')
plt.ylabel('Average Relation of Equilateral Triangles')
plt.grid('on')
plt.show()



