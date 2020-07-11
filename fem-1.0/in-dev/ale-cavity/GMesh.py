# -*- coding: utf-8 -*-
import scipy as sp

# Escrito por Luís Cunha
# Última atualização : 27/10/2017

# Código para ler aqrquivos de malha gerados pelo software Gmsh:

# $MeshFormat
# 2.2 0 8

# Code adapted from: https://github.com/jaabell/gmshtranslator

# Lista para referência:


# .dirichlet_points -- valor de dirichlet nos pontos
# .neumann_points -- valor de dirichlet nos pontos
# .Nphys -- numero de grupos físicos
# .Nnodes -- numero de pontos
# .Nelem -- número de elementos (incluindo elementos de fronteira e "meio da malha")
# .physical_groups -- indices dos grupos físicos na malha

# .nodes_in_physical_groups -- ATUALIZAR
# .physical_group_dims -- ATUALIZAR
# .physical_group_names -- ATUALIZAR

# .X , .Y, .Z -- vetores de posição das respectivas coordenadas
# .IEN -- matriz de conectividade


class GMesh:
    def __init__(self, mshfilename):

        self.mshfilename = mshfilename
        self.mshfid = open(mshfilename, "r")

        # Marcadores para ler as respectivas partes do .msh
        reading_physnames = 0
        reading_nodes = 0
        reading_elements = 0
        dirichlet_count = 0

        self.Nphys = 0
        self.Nnodes = 0
        self.Nelem = 0
        self.physical_groups = []
        self.physical_groups_dirichlet = []
        self.physical_groups_neumann = []
        self.nodes_in_physical_groups = {}
        self.physical_group_dims = {}
        self.physical_group_names = {}

        linenumber = 1
        for line in self.mshfid:

            # ----------------------------------------------------------------
            # Identificando as seções do arquivo (PhysGroup, Nodes, Elements, ...)
            # ---------------------------------------------------------------

            # Começo de cada seção
            if line.find("$PhysicalNames") >= 0:
                reading_physnames = 1
                continue
            if line.find("$Nodes") >= 0:
                reading_nodes = 1
                continue
            if line.find("$Elements") >= 0:
                reading_elements = 1
                continue

            # Fim de cada seção
            if line.find("$EndPhysicalNames") >= 0:
                reading_physnames = 0
                continue
            if line.find("$EndElements") >= 0:
                reading_elements = 0
                continue
            if line.find("$EndNodes") >= 0:
                reading_nodes = 0
                continue

            # ----------------------------------------------------------------
            # Lendo as seções do arquivo
            # ---------------------------------------------------------------

            # Número de pontos, elementos e physical groups
            if reading_physnames == 1:
                self.Nphys = sp.int32(line)
                reading_physnames = 2
                continue
            if reading_nodes == 1:
                self.Nnodes = sp.int32(line)
                reading_nodes = 2
                continue
            if reading_elements == 1:
                self.Nelem = sp.int32(line)
                reading_elements = 2
                continue

            # Dados de Physical Groups
            if reading_physnames == 2:
                sl = line.split()
                grpdim = sp.int32(sl[0])  # dimensão do grupo físico (0 = point, 1 = line, 2 = surface, 3 = volume)
                physgrp = sp.int32(sl[1])  # número de indicação do grupo
                grpname = (" ".join(sl[2:]))[1:-1]
                self.physical_group_dims[physgrp] = grpdim
                self.physical_group_names[physgrp] = grpname
                condition = grpname.split()
                if condition[0] == "dirichlet":
                    self.physical_groups_dirichlet.append(physgrp)
                if condition[0] == "neumann":
                    self.physical_groups_neumann.append(physgrp)
                continue

            # Dados dos vetores de posição
            if reading_nodes == 2:
                self.X = sp.zeros(self.Nnodes)
                self.Y = sp.zeros(self.Nnodes)
                self.Z = sp.zeros(self.Nnodes)
                reading_nodes = 3
            if reading_nodes == 3:
                sl = sp.array(line.split(), dtype=sp.float32)
                self.X[int(sl[0])-1] = sl[1]
                self.Y[int(sl[0])-1] = sl[2]
                self.Z[int(sl[0])-1] = sl[3]
                continue

            # Dados dos elementos e IEN
            if reading_elements == 2:
                self.IEN = sp.zeros((self.Nelem,3), dtype=sp.int32)
                self.Boundary = sp.zeros((self.Nelem,2), dtype=sp.int32)
                reading_elements = 3
            if reading_elements == 3:
                sl = sp.array(line.split(), dtype=sp.int32)

                eletag = sl[0]
                eletype = sl[1]
                ntags = sl[2]
                physgrp = 0

                if ntags >= 2:
                    physgrp = sl[3]
                    nodelist = sl[(3 + ntags)::]

                    if physgrp in self.physical_groups:
                        self.nodes_in_physical_groups[physgrp][nodelist] = 1
                    else:
                        self.nodes_in_physical_groups[physgrp] = -sp.ones(self.Nnodes + 1, dtype=sp.int16)
                        self.nodes_in_physical_groups[physgrp][nodelist] = 1
                        self.physical_groups.append(physgrp)
                        pass

                if eletype == 2:
                    self.IEN[sl[0]-1] = [sl[3+ntags]-1, sl[4+ntags]-1, sl[5+ntags]-1]
                if eletype == 1:
                    self.Boundary[sl[0]-1] = [sl[3+ntags]-1, sl[4+ntags]-1]

            linenumber += 1
        #Fim do loop (for line in self.mshfid:)

        # --------------------------------------------------------------------
        #       IEN and Boundary elements
        # --------------------------------------------------------------------

        listaien = []
        listabound = []
        for i in range(self.Nelem):
            if sp.all(self.IEN[i] == 0):
                listaien.append(i)

        for i in range(self.Nelem):
            if sp.all(self.Boundary[i] == 0):
                listabound.append(i)

        self.IEN = sp.delete(self.IEN, listaien, axis=0)
        self.Boundary = sp.delete(self.Boundary, listabound, axis=0)


        # --------------------------------------------------------------------
        #  Condição de DIRICHLET
        # --------------------------------------------------------------------

        dirichlet_count = 0
        for i in self.physical_groups_dirichlet:
            for j in range(1, self.Nnodes+1):
                linha = self.physical_group_names[i].split()
                if self.nodes_in_physical_groups[i][j] == 1:
                    dirichlet_count += 1

        self.dirichlet_points = sp.zeros((dirichlet_count, 2))
        self.Boundary_Nodes = sp.zeros(dirichlet_count, dtype=int)
        counter = 0
        for i in self.physical_groups_dirichlet:
            linha = self.physical_group_names[i].split()
            for j in range(1, self.Nnodes+1):
                if self.nodes_in_physical_groups[i][j] == 1:
                    self.Boundary_Nodes[counter] = int(j)-1
                    self.dirichlet_points[counter][0] = j
                    self.dirichlet_points[counter][1] = linha[1]
                    counter += 1

        # --------------------------------------------------------------------
        #  Condição de NEUMANN
        # --------------------------------------------------------------------

        boundarysize = len(self.Boundary)
        self.Boundary_Neumann = sp.zeros((boundarysize, 2), dtype=sp.int32)

        neumann_count = 0
        for i in self.physical_groups_neumann:
            for j in range(1, self.Nnodes+1):
                linha = self.physical_group_names[i].split()
                if self.nodes_in_physical_groups[i][j] == 1:
                    neumann_count += 1

        self.neumann_points = sp.zeros((neumann_count, 2))
        counter = 0
        for i in self.physical_groups_neumann:
            linha = self.physical_group_names[i].split()
            for j in range(1, self.Nnodes+1):
                if self.nodes_in_physical_groups[i][j] == 1:
                    self.neumann_points[counter][0] = j
                    self.neumann_points[counter][1] = linha[1]
                    counter += 1

        for i in range(boundarysize):
            mark1 = 0
            mark2 = 0
            for j in range(len(self.neumann_points)):
                if self.Boundary[i,0] == self.neumann_points[j,0]:
                    self.Boundary_Neumann[i, 0] = self.Boundary[i,0]
                    mark1 = 1
                if self.Boundary[i, 1] == self.neumann_points[j, 0]:
                    self.Boundary_Neumann[i, 1] = self.Boundary[i, 1]
                    mark2 = 1
            if mark1 == 0 or mark2 == 0:
                self.Boundary_Neumann[i, 0] = 0
                self.Boundary_Neumann[i, 1] = 0

        lista_neu = []
        for i in range(boundarysize):
            if sp.all(self.Boundary_Neumann[i] == 0):
                lista_neu.append(i)
        self.Boundary_Neumann = sp.delete(self.Boundary_Neumann, lista_neu, axis=0)


#----------------------------------------------------------------------
#-------------------- Other Mesh functions ----------------------------
#----------------------------------------------------------------------

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
