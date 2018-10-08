import numpy as np
from scipy.linalg import toeplitz
from scipy.sparse import csr_matrix as sparse


def buildD2(n):
    # 1D 2nd order centered finite difference operator
    # from n+1 nodes to n cell-centers
    # D : n x (n+1)
    c=np.zeros(n)
    c[0] = -1
    r=np.zeros(n+1)
    r[0:2] = [-1, 1]
    hD = toeplitz(c,r)
    return sparse(hD)
