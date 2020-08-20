# -*- coding: utf-8 -*-

import numpy as np
import meshio

#from cython_func.neighbours.cy import cy_neighbours
from neighbours import *



class Mesh:
    def __init__(self, mshfile):

        _msh = meshio.read(mshfile)

        self.x = _msh.points[:,0]
        self.y = _msh.points[:,1]
        self.ien = _msh.cells_dict['triangle']

        self.neighbour_nodes = None
        self.neighbour_elements = None

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

        self.boundary_names = [[] for i in range(self.num_nodes)]

        for elem in range(len(self.ien_boundary)):
            self.boundary_names[self.ien_boundary[elem][0]] = \
                self.ien_boundary_types[elem]
            self.boundary_names[self.ien_boundary[elem][1]] = \
                self.ien_boundary_types[elem]

        self.dirichlet_nodes = None

    def get_boundary_with_name(self, name):
        ''' Returns all boundary points that have the same type
            *arg name: Type of the boundary (ex: 'wall', 'outlet')
        '''
        points = []
        for i in range(len(self.boundary_names)):
            if self.boundary_names[i] == name:
                points.append(self.boundary_nodes[i])
        return points


    def set_dirichlet_nodes(self, _boundaryTypes):
        ''' Sets a list with all boundary nodes that have a
        prescribed dirichlet condition.
            *arg _boundaryTypes: list of boundary names ex.: ['inlet', 'wall']
        '''
        if self.dirichlet_nodes is not None:
            return
        _dirichlet_nodes = []
        for type in _boundaryTypes:
            _dirichlet_nodes.extend(self.get_boundary_with_name(type))
            print(_dirichlet_nodes)
        # self.dirichlet_nodes = list(np.unique(_dirichlet_nodes))
        # self.dirichlet_nodes = np.unique(_dirichlet_nodes)
        self.dirichlet_nodes = _dirichlet_nodes


    def set_neighbours(self, _file=None):
        ''' Sets the neighbours Data Arrays
        returns: Neighbour Elements, Neighbour Nodes
        *kwarg _file: If given, uses the file contents as the structures
                        instead of calculating
        '''
        if _file is not None:
            f = open(_file+'-neighbourNodes.txt', 'r')
            self.neighbour_nodes = []
            for line in f:
                temp = line[1:-2].split(', ')
                for i in range(len(temp)):
                    temp[i] = int(temp[i])
                self.neighbour_nodes.append(temp)
            f.close()

            f = open(_file+'-neighbourElements.txt', 'r')
            self.neighbour_elements = []
            for line in f:
                temp = line[1:-2].split(', ')
                for i in range(len(temp)):
                    temp[i] = int(temp[i])
                self.neighbour_elements.append(temp)
            f.close()
        else:
            #self.neighbour_nodes, self.neighbour_elements = cy_neighbours()
            self.neighbour_elements, self.neighbour_nodes = py_neighbours(self.num_nodes, self.ien)


#----------------------------------------------------------------------
#-------------------- Other Mesh functions ----------------------------
#----------------------------------------------------------------------


def distance(_a,_b):
    size = len(_a)
    sum = 0
    for i in range(size):
        sum += (_b[i] - _a[i])**2
    return sp.sqrt(sum)


def smoothMesh(_mesh, _dt):


    if _mesh.neighbour_nodes is None:
        _mesh.set_neighbours()

    xx = sp.copy(_mesh.x)
    yy = sp.copy(_mesh.y)
    vx_disp = sp.zeros(len(xx))
    vy_disp = sp.zeros(len(xx))
    for i in range(len(_mesh.neighbour_nodes)):

        flag = 0
        for k in _mesh.boundary_nodes:
            if i == k:
                flag = 1
                break

        if flag == 1:
            continue

        vertex_position = sp.array([xx[i], yy[i]])
        nghN = _mesh.neighbour_nodes[i]
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
