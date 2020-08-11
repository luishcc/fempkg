# -*- coding: utf-8 -*-

import numpy as np
import meshio


class Mesh:
    def __init__(self, mshfilename):

        _msh = meshio.read(mshfilename)

        self.x = _msh.points[:,0]
        self.y = _msh.points[:,1]
        self.ien = _msh.cells_dict['triangle']

        self.num_nodes = len(self.x)
        self.num_elem = len(self.ien)

        self.boundary_TypeTags = \
            list(_msh.cell_data_dict['gmsh:physical']['line']-1)

        self.boundary_Types = list(_msh.field_data.keys())

        self.ien_boundary = _msh.cells_dict['line']

        self.ien_boundary_types = \
            [self.boundary_Types[e] for e in self.boundary_TypeTags]

        self.boundary_nodes = np.unique(
            self.ien_boundary.reshape(self.ien_boundary.size) )

        self.boundary_names = [[] for i in range(len(self.boundary_nodes))]

        for elem in range(len(self.ien_boundary)):
            self.boundary_names[self.ien_boundary[elem][0]] = \
                self.ien_boundary_types[elem]
            self.boundary_names[self.ien_boundary[elem][1]] = \
                self.ien_boundary_types[elem]


    def get_boundary_with_name(self, name):
        points = []
        for i in range(len(self.boundary_names)):
            if self.boundary_names[i] == name:
                points.append(self.boundary_nodes[i])
        return points


#----------------------------------------------------------------------
#-------------------- Other Mesh functions ----------------------------
#----------------------------------------------------------------------


def distance(_a,_b):
    size = len(_a)
    sum = 0
    for i in range(size):
        sum += (_b[i] - _a[i])**2
    return sp.sqrt(sum)


def smoothMesh(_neighbour_nodes, _mesh, _x, _y, _dt):

    xx = sp.copy(_x)
    yy = sp.copy(_y)
    vx_disp = sp.zeros(len(xx))
    vy_disp = sp.zeros(len(xx))
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
            _index = nghN[j]
            distance_vectors[j][0] = xx[_index] - vertex_position[0]
            distance_vectors[j][1] = yy[_index] - vertex_position[1]
            displacement_vector += distance_vectors[j]
        displacement_vector *= (1./num_nghb)
        displacement_velocity = displacement_vector / _dt
        vx_disp[i] = displacement_velocity[0]
        vy_disp[i] = displacement_velocity[1]

      #  print nghN, "\n", distance_vectors, "\n", displacement_vector, "\n", new_position

    return vx_disp, vy_disp



def weighted_smoothMesh(_neighbour_nodes, _boundary, _x, _y, _dt):

    xx = sp.copy(_x)
    yy = sp.copy(_y)
    vx_disp = sp.zeros(len(xx))
    vy_disp = sp.zeros(len(xx))
    for i in range(len(_neighbour_nodes)):

        flag = 0
        for k in _boundary:
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
        sum_length = 0
        for j in range(num_nghb):
            _index = nghN[j]
            length = distance([xx[_index], yy[_index]], vertex_position)
            distance_vectors[j][0] = xx[_index] * length
            distance_vectors[j][1] = yy[_index] * length
            sum_length +=  length
            displacement_vector += distance_vectors[j]
        displacement_vector = displacement_vector/sum_length - vertex_position
        displacement_velocity = displacement_vector / _dt
        vx_disp[i] = displacement_velocity[0]
        vy_disp[i] = displacement_velocity[1]


    return vx_disp, vy_disp
