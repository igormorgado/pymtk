from matrizestencil import banded
from numpy import concatenate


def PiTilde(Pi, Nu, s):
    """Returns the PiTilde matrix given the Pi matrix, Nu Nullspace
    and the 's' stencil"""
    missingRowsCount = s.shape[0] - Pi.shape[0]
    A = banded(s, s.shape[0]-1)[-missingRowsCount:]
    PiT = concatenate((Pi, A, Nu)).T
    return PiT
