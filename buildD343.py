import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix as sparse
from scipy.sparse import vstack
from scipy.linalg import toeplitz


def buildD343(N):
    # construct 4th order DIV matrix : (N+1) -> N

    d343 = np.array([(-199)/208, 271/312, 7/52, (-5)/104, 1/624])
    q343 = np.array([13/12, 7/8, 25/24])
    s4 = np.array([1/24,-9/8,9/8,-1/24])

    n = len(s4)
    r4 = sparse((s4, (np.zeros(n), np.arange(n))), shape=[1, N+1])
    c4 = sparse(([s4[0]], ([0], [0])), shape=[N-2, 1])

    lnd = len(d343)

    rd343 = sparse((d343, (np.zeros(lnd), np.arange(lnd))), shape=(1,N+1))

    D = vstack([
            rd343,
            sparse(toeplitz(c4.todense(), r4.todense())),
            -rd343[::-1,::-1]])

    lnq = len(q4)
    q = np.ones(N)
    q[:lnq] = q4
    q[-lnq:] = q4[::-1]
    Q = np.diag(q)

    return D, sparse(Q)
