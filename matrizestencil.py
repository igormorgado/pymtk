# Algoritmo 3 - Cria matriz bandada do estencil MDF

from numpy import asarray, zeros, concatenate
from numpy.lib.stride_tricks import as_strided as strided
from scipy.sparse import diags

def MatrizEstencil(g, N):
    n=len(g)
    T = zeros((N,N+n-1))
    for x in range(N):
        T[x][x:x+n]=g
    return T


def MatrizEstencilEsparsa(g, N):
    n=len(g)
    return diags(g, range(N), shape=(N,N+n-1)).toarray()


def MatrizEstencilRapida(g, N):
    """Creates a `g` generated banded matrix with 'N' rows."""
    g = asarray(g)
    p = zeros(N-1, dtype=g.dtype)
    b = concatenate((p, g, p))
    s = b.strides[0]
    n = len(g)
    return strided(b[N-1:], shape=(N, n+N-1), strides=(-s, s))
