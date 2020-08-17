# -*- coding: utf-8 -*-
import scipy as sp



class Read:
    def __init__(self, mshfilename):

        self.mshfilename = mshfilename
        self.mshfid = open(mshfilename, "r")

        reading_physnames = 0
        reading_nodes = 0
 
        self.T = []
        self.V = []

        linenumber = 1
        for line in self.mshfid:

            # ----------------------------------------------------------------
            # Identificando as seções do arquivo (PhysGroup, Nodes, Elements, ...)
            # ---------------------------------------------------------------

	    sl = line.split()

            # Começo de cada seção
            if sl[0] == "T":
                reading_physnames = 1
                continue
            if line.find("V") >= 0:
                reading_nodes = 1
                continue
            

            # Fim de cada seção
            if sl[0] == "EndT":
                reading_physnames = 0
                continue

            if line.find("EndV") >= 0:
                reading_elements = 0
                continue
           

            # ----------------------------------------------------------------
            # Lendo as seções do arquivo
            # ---------------------------------------------------------------

            # Dados de T
            if reading_physnames == 1:
                sl = line.split()
                self.T.append(float(sl[0]))
                continue

            # Dados de V
            if reading_nodes == 1:
                sl = line.split()
                self.V.append(float(sl[0]))
                continue

