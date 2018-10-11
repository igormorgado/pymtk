# Algoritmo 5 - Vetor gerador para matriz de Vandermonde no
#               contorno do dominio para o operador gradiente

from numpy import linspace, zeros, array


def GeradorContornoGradiente(k):
    """Returns a matriz with 'k' generators for divergent near boundaries"""
    l = int((3/2) * k)
    g = []
    for i in range(int(l/2)):
        # These are the points near the boundary
        g_ = zeros(l)
        g_[0] = 0
        g_[1:] = linspace(1/2, (l-.5), l-1, endpoint=False)
        g.append(g_ - i)

    for i in range(int(l/2), k):
        # These are the interior points
        g.append(linspace(-1/2, (l-.5), l, endpoint=False)-(i-1))

    return array(g, dtype=float)
