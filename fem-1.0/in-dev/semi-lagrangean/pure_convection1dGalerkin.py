# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import libmalha as mlh

# Escrito por Luís Cunha
# Última atualização : 19/03/2019


# Código para a solução da EDP DT
#                              --- = 0
#                              Dt


#

L = 1                       # Comprimento da Malha 1D
num_ele = 100                  # Número de elementos
tempo = 100000
dt = 0.00001
v = 0.1
n_points = num_ele+1
pi10 = np.pi*10

    # Definição da malha

x, ien = mlh.lin1d(num_ele, L)

dx = np.zeros(num_ele)
dxa = 0
for i in range(num_ele):
    dx[i]= x[i+1]-x[i]
    dxa += dx[i]

dxa = dxa/float(num_ele)

     # Matrizes de massa, rigidez e carregamento por elemento

m = np.array([[1,2], [2,1]])
g = np.array([[-0.5,0.5], [-0.5,0.5]])


    # Montagem das Matrizes globais

M = np.zeros((num_ele+1, num_ele+1))
G = np.zeros((num_ele+1, num_ele+1))

for elem in range(num_ele):
    for i_local in range(2):
        i_global = ien[elem, i_local]-1
        for j_local in range(2):
            j_global = ien[elem, j_local] -1
            M[i_global, j_global] = M[i_global, j_global] + (dxa / 6.) * m[i_local, j_local]
            G[i_global, j_global] = G[i_global, j_global] + g[i_local, j_local]



LHS = M + v * G * dt


    #Implementando Condicao inicial


T_ini = np.zeros(n_points)

for i in range(n_points):
    if x[i] <= 0.1:
        T_ini[i] = np.sin(pi10*x[i])



T_i = np.copy(T_ini)

    # Solucao do sistema

interval = np.linspace(0, tempo, 101)

cont = 1
plt.figure(0)
plt.plot(x, T_ini, 'r.')
plt.xlim(-0.01, L + 0.01)
plt.ylim(-0.3, 1.3)
plt.title('Galerkin, t = 0')
plt.xlabel('x')
plt.ylabel('c')
plt.savefig('/home/luis/Pictures/imagem2/fig' + str(0))

# print interval
# exit()


for i in range(0, tempo):
    B = np.dot(M, T_i)
    T_d = np.linalg.solve(LHS, B)
    T_i = np.copy(T_d)

    print i*0.001, ' %'
    step = dt * i

    if i == int(interval[cont]):
        plt.figure(i+1)
        plt.plot(x, T_d, 'r.')
        plt.xlim(-0.01, L+0.01)
        plt.ylim(-0.3, 1.3)
        plt.title('Galerkin, t = '+str(step))
        plt.xlabel('x')
        plt.ylabel('c')
        plt.savefig('/home/luis/Pictures/imagem2/fig'+str(i+1))
        cont += 1
        pass
    else:
        continue
