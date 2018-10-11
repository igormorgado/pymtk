import numpy as np

def VetorSeletor(N, n=1):
    """Returns the N-dimensional selector vector to the n-th derivative"""
    o = np.zeros(N)
    o[n] = 1
    return o 
