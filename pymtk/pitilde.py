# Algoritmo 7 - Constroe a matriz \tilde(\Pi)

from matrizestencil import MatrizEstencil
from numpy import concatenate


def MatrizPiTilde(Pi, Nu, s):
    """Returns the PiTilde matrix given the Pi matrix, Nu Nullspace
    and the 's' stencil"""
    missingRowsCount = s.shape[0] - Pi.shape[0]
    N = Pi.shape[1] - s.shape[0] + 1
    A = MatrizEstencil(s, N)[-missingRowsCount:]
    PiT = concatenate((Pi, A, Nu)).T
    return PiT
