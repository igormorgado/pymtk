# Algoritmo 4 - Vetor gerador para matriz de Vandermonde no 
#               contorno do dominio para o operador divergente

from numpy import linspace


def divergent_generator(k,h):
    """Returns a matriz with 'k' generators for 
    divergent near boundaries"""
    l = int((3/2) * k)
    g0 = linspace(-h/2, h*(l-.5), l, endpoint=False)
    sub = linspace(0, k*h, k, endpoint=False)
    return g0 - sub.reshape(sub.shape[0],1)

