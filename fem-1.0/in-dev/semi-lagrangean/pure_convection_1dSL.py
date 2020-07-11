# -*- coding: utf-8 -*-
import scipy as sp
import matplotlib.pyplot as plt
import libmalha as mlh


# Escrito por Luís Cunha
# Última atualização : 19/03/2019


# Código para a solução da EDP DT
#                              --- = 0
#                              Dt

# através do semi-lagrangeano
# Serão utilizadas condições de contorno  T'(0,t) = T'(L,t) = 0 e T(x, 0) = sin()
#

L = 1                       # Comprimento da Malha 1D
num_ele = 100                  # Número de elementos
tempo = 100
dt = 0.01
n_points = num_ele+1
vx = 0.1
pi10 = 10 * sp.pi

    # Definição da malha

x_i, ien = mlh.lin1d(num_ele, L)
ien -= 1
x_d = x_i - vx * dt


    #Implementando Condicao inicial

T_ini = sp.zeros(n_points)

for i in range(n_points):
    if x_i[i] <= 0.1:
        T_ini[i] = sp.sin(pi10*x_i[i])


def Linear1D_v2(_npoints, _nelem, _IEN, _xn, _vx, _dt, _scalar):
    xd = _xn - _vx * _dt

    scalar = sp.zeros([_npoints, 1], dtype=float)

    for i in range(0, _npoints):
        x = float(xd[i])

        breaking = 0
        length = []

        for e in range(0, _nelem):
            v1 = _IEN[e][0]
            v2 = _IEN[e][1]

            x1 = float(_xn[v1])
            x2 = float(_xn[v2])

            len1 = abs(x2 - x)
            len2 = abs(x1 - x)
            lent = abs(x1 - x2)

            Li = len1 / lent
            Lj = len2 / lent

            alpha = [Li, Lj]
            alpha = sp.array(alpha)

            if sp.all(alpha >= 0.0) and sp.all(alpha <= 1.0):
                Ni = Li
                Nj = Lj

                scalar1 = _scalar[v1]
                scalar2 = _scalar[v2]

                scalar[i] = Ni * scalar1 + Nj * scalar2
                breaking = 1
                break

            else:
                x_a = x1 - x
                x_b = x2 - x

                length1 = sp.sqrt(x_a ** 2)
                length2 = sp.sqrt(x_b ** 2)

                a_1 = [v1, length1]
                a_2 = [v2, length2]

                length.append(a_1)
                length.append(a_2)

        if breaking == 0:
            length_min = min(length, key=lambda k: k[1])
            node = length_min[0]
            scalar[i] = _scalar[node]

    return scalar

    #Calculando


T_i = sp.copy(T_ini)


for i in range(tempo):

    T_d = Linear1D_v2(n_points, num_ele, ien, x_i, vx, dt, T_i)
    T_i = sp.copy(T_d)

    step = dt*i
    plt.figure(i)
    plt.plot(x_i, T_i, 'r.')
    plt.xlim(-0.01, L+0.01)
    plt.ylim(-0.3, 1.3)
    plt.title('Semi-Lagrangean, t = '+str(step))
    plt.xlabel('x')
    plt.ylabel('c')
    plt.savefig('/home/luis/Pictures/imagem/fig'+str(i))