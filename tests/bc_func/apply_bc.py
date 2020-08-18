import numpy as np


def apply_bc_dirichlet(_LHS, _mesh, _bc_value)

    A = np.copy(_LHS)
    num_bc = len(_mesh.boundary_nodes)

    for i in range(num_bc):
        index = int(Boundary[i])
        value = Wcc[index]
        for j in range(nodes):
            ccomega[j] -= value * LHS[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1



"""
 def gaussian_elimination(_self, _LHS0, _neighbors_nodes):
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float)
  _self.neighbors_nodes = _neighbors_nodes

  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0

   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0

"""
