from sympy import nsimplify, Matrix


def rationalize(x, tolerance=10e-9):
    """Returns a Matrix with symbolic rational numbers"""
    try:
        Mx = Matrix(x)
    except TypeError:
        raise("Input isn't a valid matrix type. Need to be a sequence")

    shape = Mx.shape
    rationalized = [nsimplify(e, tolerance=tolerance) for e in Mx]
    return Matrix(rationalized).reshape(*shape)


