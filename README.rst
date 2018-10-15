=====
pymtk
=====


.. image:: https://img.shields.io/pypi/v/pymtk.svg
        :target: https://pypi.python.org/pypi/pymtk

.. image:: https://img.shields.io/travis/igormorgado/pymtk.svg
        :target: https://travis-ci.org/igormorgado/pymtk

.. image:: https://readthedocs.org/projects/pymtk/badge/?version=latest
        :target: https://pymtk.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Mimetic Methods for Python.

This package try to build Castillo Grone mimetic operators to apply to finite
differences models.

* Free software: MIT license
* Documentation: https://pymtk.readthedocs.io.


How to use
----------

.. code-block:: python

import pymtk
import numpy as np

# Number or discretized points
N = 11

# Operator order
k = 4

# Fourth order Divergent and Gradient
D_4 = pymtk.Divergent(order=k)
G_4 = pymtk.Gradient(order=k)

# Weight Matrices
Q = np.diag(D_4.weight_vector(N))
P = np.diag(G_4.weight_vector(N))

# Then you can call the operator passing the number of
# grid points and it will return a numpy matrix
D = D_4(N) 

::



Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
