# Algoritmo 4 - Vetor gerador para matriz de Vandermonde no 
#               contorno do dominio para o operador divergente

from numpy import linspace


def GeradorContornoDivergente(k):
    """Returns a matriz with 'k' generators for 
    divergent near boundaries"""
    l = int((3/2) * k)
    g0 = linspace(-1/2, (l-.5), l, endpoint=False)
    sub = linspace(0, k, k, endpoint=False)
    return g0 - sub.reshape(sub.shape[0],1)

