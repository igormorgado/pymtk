# Algoritmo 8 - Calcula vetor mimetico da submatriz A(alpha)

from numpy import zeros, float_


def VetorMimetico(s):
    """Returns the resultant mimetic vector for the submatrix A"""

    k = len(s)
    l = 3*k//2

    h = zeros(l,dtype=float_)
    h[0] = 1

    i = k//2 + 1
    h[i:]=s.cumsum()[:(l-i)]

    return -h
