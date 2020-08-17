import cython
cimport cython
import numpy as np
cimport numpy as np

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_global = np.zeros((_numnode, _numnode), dtype="float64")
    m_global = np.zeros((_numnode, _numnode), dtype="float64")
    gx_global = np.zeros((_numnode, _numnode), dtype="float64")
    gy_global = np.zeros((_numnode, _numnode), dtype="float64")

    _matrix(_x, _y, _numele, _numnode, _ien,
          k_global, m_global,
          gx_global,  gy_global)

    return  k_global, m_global, gx_global, gy_global


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _matrix(double[:] _x, double[:] _y, int _numele, int _numnode,
        unsigned long[:,:] _ien,
        double[:,:] k_global, double[:,:] m_global,
        double[:,:] gx_global,  double[:,:] gy_global):

    cdef:
        unsigned long elem, i, j, i_global, j_global
        int i_local = 0, j_local = 0
        double Area
        double[3] a, b, c, xx, yy
        double[3][3] k_local, gx_local, gy_local
        double[3][3] m_local = [[2, 1, 1], [1, 2, 1], [1, 1, 2]]

    for elem in range(_numele):
        i = 0
        for i in range(3):
            idx = _ien[elem][i]
            xx[i] = _x[idx]
            yy[i] = _y[idx]

        a[0] = xx[0] * yy[1] - xx[1] * yy[0]
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        Area = (a[0] + a[1] + a[2]) / 2.

        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]

        i = 0
        for i in range(3):
            j = 0
            for j in range(3):
                k_local[i][j] = (b[i] * b[j] + c[i] * c[j]) / (4 * Area)
                gx_local[i][j] = b[j] * (1/6.)
                gy_local[i][j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem][i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global][j_global] += k_local[i_local][j_local]
                m_global[i_global][j_global] += m_local[i_local][j_local] *\
                                                (Area/12.)
                gx_global[i_global][j_global] += gx_local[i_local][j_local]
                gy_global[i_global][j_global] += gy_local[i_local][j_local]
