# -*- coding: utf-8 -*-

"""Main module."""

import abc

from numpy import ( arange, vander, zeros, array, newaxis, concatenate, float_,
                  linspace, flipud, fliplr, ones )
from sympy import linsolve, Matrix
from scipy.linalg import lstsq


class MimeticOperator:
    
    def __init__(self, order):
        assert (order % 2 == 0 and order > 0), "Order must be an even number > 0"

        self.order = order
        self.l = 3*self.order//2

        #self._N = None

        # Build the operator
        self._stencil()             # s
        self._mimetic_vector()      # d
        self._boundary_stencil()    # Pi, Nu
        self._pitilde()             # PiTilde
        self._weights()             # Q and lambda
        self._boundary_rows()       # A(alpha)


    @staticmethod
    def sliding_window(window, elements):
        n=len(window)
        T = zeros((elements,elements+n-1))
        for x in range(elements):
            T[x][x:x+n]=window
        return T

    @staticmethod
    def order_selector_vector(N, n=1):
        """Returns the N-dimensional selector vector to the n-th derivative"""
        o = zeros(N)
        o[n] = 1
        return o


    def _generator_vector(self):
        return arange((1-self.order)/2, self.order/2, 1)

    @property
    def generator_vector(self):
        return self._generator_vector()


    def _stencil(self):
        g = self.generator_vector
        T = vander(g, increasing=True).T
        o = self.order_selector_vector(len(g), 1)
        s = lstsq(T, o)[0]
        #import ipdb; ipdb.set_trace()
        self.stencil =  s


    @abc.abstractmethod
    def _boundary_stencil(self):
        raise NotImplementedError('operator must implement boundary stencil')


    def _pitilde(self):
        """Returns the PiTilde matrix given the Pi matrix, Nu Nullspace
        and the 's' stencil"""
        missingRowsCount = self.order - self.Pi.shape[0]
        n = self.Pi.shape[1] - self.order + 1
        A = self.sliding_window(self.stencil, n)[-missingRowsCount:]
        PiT = concatenate((self.Pi, A, self.Nu)).T
        self.pitilde = PiT


    def _mimetic_vector(self):
        h = zeros(self.l, dtype=float_)
        h[0] = 1
        i = self.order//2 + 1
        h[i:] = self.stencil.cumsum()[:(self.l-i)]

        self.mimetic_vector = -h


    def _weights(self):
        solution = lstsq(self.pitilde, self.mimetic_vector)[0]
        i = self.order//2 - 1
        self.weights = solution[:-i]
        self.lambda_ = solution[-i:]

    
    def _boundary_rows(self):
        """Returns operator boundary rows (aka submatrix A)"""
        A = self.pitilde.T[:self.order]

        alpha = self.lambda_ / self.weights[:len(self.lambda_)]
        for i in range(len(alpha)):
            A[i] += alpha[i] * self.Nu[i]

        self.boundary_rows = A


    def _weight_vector(self, N):
        W = ones(N)
        W[:self.order] = self.weights
        W[-self.order:] = self.weights[::-1]
        return W


    def weight_vector(self, N):
        return self._weight_vector(N)


    def _operator_matrix(self):
        """To build the operator matrix we need the number of discretized points"""


    @abc.abstractmethod
    def __call__(self, N=None):
        raise NotImplementedError('operator must implement boundary stencil')


    

class Divergent(MimeticOperator):
    """Creates the Stencil near the boundaries

    Where:
        k: Order of accuracy

    Returns: PI, NU
        Mimetic submatrix at boundary
        Kernel of boundary Vandermonde matrix
    """

    def __init__(self, order):
        self.order = order
        super().__init__(self.order)

    def boundary_generator(self):
        """Returns a matriz with 'k' generators for 
        divergent near boundaries"""
        g = linspace(-.5, (self.l-.5), self.l, endpoint=False)
        s = linspace(0, self.order, self.order, endpoint=False)
        return g - s[:, newaxis]

    def _boundary_stencil(self):
        o = self.order_selector_vector(self.order+1)
        d = self.order//2 - 1
        g = self.boundary_generator()
        Pi = zeros([d, self.l])

        lNu = []
        for i in range(d):
            V_tmp = vander(g[i], self.order+1, increasing=True).T
            r_tmp = linsolve(Matrix(concatenate((V_tmp, o[:,newaxis]), axis=1)))
            r_tmp = r_tmp.subs([(s, 0) for s in r_tmp.free_symbols])
            Pi[i] = array(*[e for e in r_tmp]).astype(float_)
            lNu.append(Matrix(V_tmp).nullspace())

        Nu = zeros([len(lNu[0]), self.l])
        for i, n in enumerate(lNu[0]):
            Nu[i] = array(n.T).astype(float_)[0]

        self.Pi = Pi
        self.Nu = Nu


    def __call__(self, N=None):
        asrmsg = f"Discretized grid must be larger than {3*self.order - 1}"
        assert (N >= 3*self.order -1), asrmsg

        # Works for Divergent N+1 -> N  matrix
        k, l = self.boundary_rows.shape

        O = zeros((N,N+1))

        # THe N-2 here ONLY works to 4th order.
        O[1:-1] = self.sliding_window(self.stencil, N-2)
        O[:k,:l] = self.boundary_rows.copy()
        O[-k:,-l:] = - flipud(fliplr(self.boundary_rows))
        return O


class Gradient(MimeticOperator):
    """Returns the Stencil near the boundaries to Gradient

    Where:
        k: Order of accuracy

    Returns: PI, NU
        Mimetic submatrix at boundary
        Kernel of boundary Vandermonde matrix
    """

    def __init__(self, order):
        self.order = order
        super().__init__(self.order)

    def boundary_generator(self):
        """Returns a matriz with 'k' generators for 
        divergent near boundaries"""
        g = []
        for i in range(self.l//2):
            # These are the points near the boundary
            g_ = zeros(self.l)
            g_[0] = 0
            g_[1:] = linspace(.5, (self.l-.5), self.l-1, endpoint=False)
            g.append(g_ - i)

        for i in range(self.l//2, self.order):
            # These are the interior points
            g.append(linspace(-1/2, (self.l-.5), self.l, endpoint=False)-(i-1))

        return array(g, dtype=float_)


    def _boundary_stencil(self):
        o = self.order_selector_vector(self.order+1)
        d = self.order//2
        g = self.boundary_generator()

        Pi = zeros([d, self.l])

        lNu = []
        for i in range(d):
            V_tmp = vander(g[i], self.order+1, increasing=True).T
            r_tmp = linsolve(Matrix(concatenate((V_tmp, o[:,newaxis]), axis=1)))
            r_tmp = r_tmp.subs([(s, 0) for s in r_tmp.free_symbols])
            Pi[i] = array(*[e for e in r_tmp]).astype(float_)
            lNu.append(Matrix(V_tmp).nullspace())

        # Aqui tem que retornar a lista de array com o nullspace completo
        Nu = zeros([len(lNu[0]), self.l])
        for i, n in enumerate(lNu[0]):
            Nu[i] = array(n.T).astype(float_)[0]

        self.Pi = Pi
        self.Nu = Nu


    def weight_vector(self, N):
        return self._weight_vector(N+1)

        
    def __call__(self, N=None):
        """ Maps N+2 --> N+1"""

        asrmsg = f"Discretized grid must be larger than {3*self.order - 1}"
        assert (N >= 3*self.order -1), asrmsg

        # Works for Divergent N+1 -> N  matrix
        k, l = self.boundary_rows.shape

        O = zeros((N+1,N+2))

        # THe N-1 here ONLY works to 4th order.
        O[1:-1] = self.sliding_window(self.stencil, N-1)
        O[:k,:l] = self.boundary_rows.copy()
        O[-k:,-l:] = - flipud(fliplr(self.boundary_rows))
        return O

