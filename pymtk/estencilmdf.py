# Algoritmo 2 - Calcula o estencil do MDF pelas matrizes de Vandermonde

from numpy import vander
from scipy.linalg import lstsq
from vetorseletor import VetorSeletor


def EstencilMDF(g):
    # Matriz de vandermonde associada
    T = vander(g, increasing=True).T

    # Vetor seletor, para primeira derivada
    o = VetorSeletor(len(g))

    # Estencil de diferencas finitas centrais da primeira derivada
    s = lstsq(T, o)[0]

    return s

