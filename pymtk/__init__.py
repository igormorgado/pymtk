# -*- coding: utf-8 -*-

"""Top-level package for pymtk."""

from sympy import nsimplify, Matrix

from pymtk.pymtk import MimeticOperator, Gradient, Divergent
import pymtk.ref as ref

author__ = """Igor Morgado"""
email__ = 'morgado.igor@gmail.com'
version__ = '0.1.0'


def rationalize(x, tolerance=10e-9):
    """Returns a Matrix with symbolic rational numbers"""
    try:
        Mx = Matrix(x)
    except TypeError:
        raise("Input isn't a valid matrix type. Need to be a sequence")

    shape = Mx.shape
    rationalized = [nsimplify(e, tolerance=tolerance) for e in Mx]
    return Matrix(rationalized).reshape(*shape)

def flipall(a):
    """Flip all axes of a given matrix (sparse or not)"""
    sl = slice(None, None, -1)
    return a[tuple([sl]*a.ndim)]
