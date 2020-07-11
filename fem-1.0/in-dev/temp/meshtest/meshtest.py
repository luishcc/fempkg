import GMesh as gm

arquivo = "malhateste"

print "Reading .msh file"
malha = gm.GMesh(arquivo+".msh")
convec = gm.GMesh(arquivo+"-conv.msh")


print malha.neumann_points
print malha.Boundary_Neumann
print malha.neumann_element_physgrp
