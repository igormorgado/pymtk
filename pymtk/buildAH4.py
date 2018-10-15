import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix as sparse
from scipy.linalg import toeplitz


def buildAH4(M):
    # the matrix A_H from Yefet and Petropoulos, ’A Staggered Fourth-Order
    # Accurate Explicit Finiite Difference Scheme for the Time-Domain Maxwell’s Equations’,
    # Computational Physics, 168(2001), pp286-315

    N = M - 1

    d4 = np.array([-23/24, 21/24, 3/24 ,-1/24])
    s4 = np.array([1/24,-9/8,9/8,-1/24])
    n = len(s4)

    r4 = sparse((s4, (np.zeros(n), np.arange(n))), shape=[1, N+1])
    c4 = sparse(([s4[0]], ([0], [0])), shape=[N-2, 1])
    lnd = len(d4)
    rd4 = sparse((d4, (np.zeros(lnd), np.arange(lnd))), shape=[1, N+1])
    D = sparse(np.concatenate((rd4.todense(), toeplitz(c4.todense(),r4.todense()), np.fliplr(-rd4.todense()))))

    return D