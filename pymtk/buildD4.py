import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix as sparse
from scipy.sparse import vstack
from scipy.linalg import toeplitz


def buildD4(N):
    # construct 4th order DIV matrix : (N+1) -> N

    d4 = np.array([[(-4751)/5192, 909/1298, 6091/15576, (-1165)/5192, 129/2596, (-25)/15576]])
    q4 = np.array([649/576, 143/192, 75/64, 551/576])

    s4 = np.array([1/24,-9/8,9/8,-1/24])
    n = len(s4)

    r4 = sparse((s4, (np.zeros(n), np.arange(n))), shape=[1, N+1])
    c4 = sparse(([s4[0]], ([0], [0])), shape=[N-2, 1])

    D = np.concatenate((
            np.zeros((1,N+1)),
            toeplitz(c4.todense(), r4.todense()),
            np.zeros((1,N+1))))

    dn, dm = d4.shape
    D[:dn,:dm] = d4
    D[-dn:,-dm:] = -np.fliplr(np.flipud(d4))

    lnq = len(q4)
    q = np.ones(N)
    q[:lnq] = q4
    q[-lnq:] = q4[::-1]
    Q = np.diag(q)

    return sparse(D), sparse(Q)
