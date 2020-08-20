import MESH as mm
import neighbours as nh
from timeit import default_timer as timer


mesh = mm.Mesh('vivTest.msh')


t1 = timer()
a1, b1 = nh.py_neighbours(mesh.num_nodes, mesh.ien)
t2 = timer()
a2, b2 = nh.py_neighbours2(mesh.num_nodes, mesh.ien)
t3 = timer()

print(t2-t1)
print(t3-t2)
