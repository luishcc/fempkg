import cython
cimport cython
import numpy as np
cimport numpy as np

def neighbours(_mesh):
    result_node = np.zeros((_mesh.num_nodes,20)), dtype="int32")-1
    result_ele = np.zeros((_mesh.num_nodes,20)), dtype="int32")-1


    _cy_neighbour(result_node, result_ele, _mesh.ien, _mesh.num_elem)

    return  result_ele, result_node


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _cy_neighbour(int[:,:] _rn, int[:,:] _re, int[:,:] _ien, int _NE):

    cdef:
        int e, i, j, v1, v2, v3

    for e in range(_NE):
        v1 = _ien[e][0]
        v2 = _ien[e][1]
        v3 = _ien[e][2]
