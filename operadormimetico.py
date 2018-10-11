# Algoritmo 10 - Constroi o operador mimetico de Castillo Grone

from matrizestencil import MatrizEstencil
from numpy import fliplr, flipud

def OCG(A, s, N):
    """Given A submatrix, build the operator matrix for N points"""
    k, l = A.shape
    O = MatrizEstencil(s, N)
    O[:k,:l] = A.copy()
    O[-k:,-l:] = -flipud(fliplr(A))
    return O
    rationalize(O, 10e-2)
    rationalize(A, 10e-2)



