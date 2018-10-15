
import numpy as np
from scipy.linalg import toeplitz

def buildD2(N):
    # N x N+1  -> 
    c=np.zeros(n)
    c[0] = -1
    r=np.zeros(n+1)
    r[0:2] = [-1, 1]
    hD = toeplitz(c,r)
    return hD

def buildG2(N):
    # N+1 x N+2  -> 
    c=np.zeros(n)
    c[0] = -1
    r=np.zeros(n+1)
    r[0:2] = [-1, 1]
    hD = toeplitz(c,r)
    return hD


