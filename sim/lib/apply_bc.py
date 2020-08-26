import numpy as np


def apply_bc_dirichlet(_LHS, _mesh, _dirichlet_nodes, _bc_value):

    A = np.copy(_LHS)
    bc_RHS = np.zeros(_mesh.num_nodes)
    num_bc = len(_dirichlet_nodes)

    for i in range(num_bc):
        index = _dirichlet_nodes[i]
        value = _bc_value[index]
        for j in _mesh.neighbour_nodes[index]:
            bc_RHS[j] -= value * _LHS[j, index]
            A[index, j] = 0
            A[j, index] = 0
            A[index, index] = 1

    return A, bc_RHS
