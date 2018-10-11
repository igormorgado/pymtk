# Algoritmo 1 - Vetor Gerador dos pontos da malha intercalada

from numpy import arange


def VetorGerador(k):
    assert (k % 2 == 0 and k > 0), "Order must be an even number > 0"
    return arange((1-k)/2, k/2, 1)
