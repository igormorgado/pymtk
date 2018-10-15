# Algoritmo 6 - Calcula o estencil e o nucleo dos pontos da
#               borda do operador gradiente

from numpy import zeros, vander, array, concatenate, newaxis
from numpy import float_

from sympy import linsolve, Matrix

from vetorseletor import VetorSeletor
from geradorgradiente import GeradorContornoGradiente


def ContornoGradiente(k):
    """Returns the Stencil near the boundaries to Gradient

    Where:
        k: Order of accuracy

    Returns: PI, NU
        Mimetic submatrix at boundary
        Kernel of boundary Vandermonde matrix
    """
    # Vetor seletor de ordem da primeira derivada
    o = VetorSeletor(k+1)

    # Kernel rank of Vandermonde matrix at boundary
    d = int(k/2)

    # Constroi os pontos do operador divergente proximo ao contorno
    g = GeradorContornoGradiente(k)

    l = 3*k//2
    Pi = zeros([d, l])
    N = []
    for i in range(d):
        V_tmp = vander(g[i], k+1, increasing=True).T
        #r_tmp = linsolve(Matrix(c_[V_tmp, o]))
        r_tmp = linsolve(Matrix(concatenate((V_tmp, o[:,newaxis]), axis=1)))
        r_tmp = r_tmp.subs([(s, 0) for s in r_tmp.free_symbols])
        Pi[i] = array(*[e for e in r_tmp]).astype(float_)
        N.append(Matrix(V_tmp).nullspace())

    Nu = zeros([len(N[0]), l])
    for i,n in enumerate(N[0]):
        Nu[i] = array(n.T).astype(float_)[0]

    return Pi, Nu
