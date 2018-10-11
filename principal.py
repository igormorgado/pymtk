#!/usr/bin/env python3
################################################################
# Algoritmo X - Codigo Principal
################################################################

import sys
from sympy import pprint

from utils import rationalize
from vetorgerador import VetorGerador
from estencilmdf import EstencilMDF
from vetorseletor import VetorSeletor
from bordagradiente import ContornoGradiente
from bordadivergente import ContornoDivergente
from pitilde import MatrizPiTilde
from vetormimetico import VetorMimetico
from pesos import VetorPeso
from submatrizA import MatrizA
from operadormimetico import OCG

k = 4               # Ordem de precisao do operador
h = 1               # Espacamento da celula
N = 3*k-1           # Menor numero de pontos possivel
l = int(3/2*k)

# Obtendo os pontos vizinhos para aproximacao de Taylor
g = VetorGerador(k)

# Obtendo o estencil do MDF
s = EstencilMDF(g)

# Constroi o vetor resultado `d` (que independe do operador)
d = VetorMimetico(s)


# Obtem operador Divergente
Pi_D, Nu_D = ContornoDivergente(k)
PiT_D = MatrizPiTilde(Pi_D, Nu_D, s)
weights_D, lambda_D = VetorPeso(PiT_D, d)
D_A = MatrizA(PiT_D, Nu_D, weights_D, lambda_D)

# Obtem operador Gradiente
Pi_G, Nu_G = ContornoGradiente(k)
PiT_G = MatrizPiTilde(Pi_G, Nu_G, s)
weights_G, lambda_G = VetorPeso(PiT_G, d)
G_A = MatrizA(PiT_G, Nu_G, weights_G, lambda_G)

# Constroi a matriz dos operadores
D = OCG(D_A, s, N)
G = OCG(G_A, s, N+3)
# L = np.dot(D, G)



sys.exit(0)

# Debug info
print("\nPi_D = ")
pprint(rationalize(Pi_D))

print(f"\nNu_D = ")
pprint(rationalize(Nu_D))

print("\nPi_G = ")
pprint(rationalize(Pi_G))

print(f"\nNu_G = ")
pprint(rationalize(Nu_G))

print(f"\nd^T = ")
pprint(rationalize(d))

