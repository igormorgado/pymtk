# Resolve o sistema linear \tilde{Pi} \tilde{q} = d
from scipy.linalg import lstsq


def VetorPeso(Pi, d):
    """Returns the weights for the diagonal matrix and the lambda"""

    k = len(d)*2//3
    i = k//2-1

    solution = lstsq(Pi, d)[0]

    weights = solution[:-i]
    lambda_ = solution[-i:]
    return weights, lambda_

