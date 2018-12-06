import pymtk
import numpy as np
# pylint: disable=E501


def stencil_2nd_order():
    return np.array([-1, 1])


def stencil_4th_order():
    return np.array([1/24, -27/24, 27/24, -1/24])


def divergent_4th_order_rojas2008(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Rojas 2008, Modelling of rupture propagation using high-order mimetic finite differences

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N x N+1
    """
    D = np.zeros((N, N+1))
    D[1:N-1] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-2)
    D[0][:6] = np.array([-1045/1142, 909/1298, 201/514, -1165/5192, 129/2596, -25/15576])
    D[N-1] = -(D[0][::-1])
    return D


def gradient_4th_order_rojas2008(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Rojas 2008, Modelling of rupture propagation using high-order mimetic finite differences

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N+1 x N+2
    """
    G = np.zeros((N+1, N+2))
    G[1:N] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-1)
    G[0][:6] = np.array([-1775/528, 1790/407, -2107/1415, 1496/2707, -272/2655, 25/9768])
    G[1][:5] = np.array([16/105, -31/24, 29/24, -3/40, 1/168])
    G[-2:] = - np.fliplr(np.flipud(G[:2]))
    return G


def divergent_4th_order_rojas2013(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Rojas 2013, Low dispersive modeling of Rayleigh waves on partly staggered grids

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N x N+1
    """
    D = np.zeros((N, N+1))
    D[1:N-1] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-2)
    D[0][:6] = np.array([-4751/5192, 909/1298, 6091/15576, -1165/5192, 129/2596, -25/15576])
    D[N-1] = -(D[0][::-1])
    return D


def gradient_4th_order_rojas2013(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Rojas 2013, Low dispersive modeling of Rayleigh waves on partly staggered grids

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N+1 x N+2
    """
    G = np.zeros((N+1, N+2))
    G[1:N] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-1)
    G[0][:6] = np.array([-47888/14245, 1790/407, -14565/9768, 8997/16280, -2335/22792, 25/9768])
    G[1][:5] = np.array([16/105, -31/24, 29/24, -3/40, 1/168])
    G[-2:] = - np.fliplr(np.flipud(G[:2]))
    return G


def divergent_24_order_rojas2013(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Rojas 2013, Low dispersive modeling of Rayleigh waves on partly staggered grids

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N x N+1
    """
    D = np.zeros((N, N+1))
    D[1:N-1] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-2)
    D[0][:2] = np.array([-1, 1])
    D[1][:4] = np.array([1/23, -26/23, 26/23, -1/23])
    D[-2:] = - np.fliplr(np.flipud(D[:2]))
    return D


def gradient_24_order_rojas2013(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Rojas 2013, Low dispersive modeling of Rayleigh waves on partly staggered grids

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N+1 x N+2
    """
    G = np.zeros((N+1, N+2))
    G[1:N] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-1)
    G[0][:3] = np.array([-8/3, 3, -1/3])
    G[1][:5] = np.array([4/39, -31/26, 44/39, -1/26])
    G[-2:] = - np.fliplr(np.flipud(G[:2]))
    return G


def divergent_4th_order_castillo2000(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Castillo 2000, Fourth and sixth order conservative finite difference
        approximations of the divergence and gradient

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N x N+1
    """
    D = np.zeros((N, N+1))
    D[1:N-1] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-2)
    D[0][:6] = np.array([-4751/5192, 909/1298, 6091/15576, -1165/5192, 129/2596, -25/15576])
    D[N-1] = -(D[0][::-1])
    return D


def gradient_4th_order_castillo2000(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Castillo 2000, Fourth and sixth order conservative finite difference
        approximations of the divergence and gradient

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N+1 x N+2
    """
    G = np.zeros((N+1, N+2))
    G[1:N] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-1)
    G[0][:6] = np.array([-1152/407, 10063/3256, 2483/9768, -3309/3256, 2099/3256, 2099/3256, -697/4884])
    G[1][:6] = np.array([0, -11/12, 17/24, 3/8, -5/24, 1/24])
    G[-2:] = - np.fliplr(np.flipud(G[:2]))
    return G


def divergent_34_order_runyan2011(N):
    """Returns a 4th order castillo grone divergent operator as defined in
    Runyan2011, Novel Higher Order FInite Difference Timings Concerning the TIme dependente Maxwell Equations

    Input:
        N: Number intervals (vector dimension)
    Returns:
        Array: N x N+1
    """
    D = np.zeros((N, N+1))
    D[1:N-1] = pymtk.MimeticOperator.sliding_window(stencil_4th_order(), N-2)
    D[0][:5] = np.array([-199/208, 271/312, 7/52, -5/104, 1/624])
    D[N-1] = -(D[0][::-1])
    return D
