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

   import numpy as np
   import pymtk


   # Keep things nicer
   np.set_printoptions(suppress=True)
   np.set_printoptions(precision=4)
   np.set_printoptions(linewidth=150)
   
   # Operator Order
   k = 4

   # Discretized points
   N = 4*k-2

   # Build the 4th order mimetic divergent to 16 points grid
   Dk = pymtk.Divergent(order=k)
   D = Dk(N)
   Q = np.diag(Dk.weight_vector(N))

   # Build the 4th order mimetic gradient to 16 points grid
   Gk = pymtk.Gradient(order=k)
   G = Gk(N)
   P = np.diag(Gk.weight_vector(N))
   
   # TESTS
   tol=1e-14


   # Divergent 
   # First Mimetic Condition
   np.all((Q @ D).sum(axis=1) < tol)

   # Second Mimetic
   tfc = np.zeros(N+1)
   tfc[0], tfc[-1] = -1, 1
   np.all(((Q @ D).sum(axis=0) - tfc) < tol)
   

   # Gradient

   # Fist mimetic condition
   np.all((P @ G).sum(axis=1) < tol)

   # Second Mimetic Condition
   tfc = np.zeros(N+2)
   tfc[0], tfc[-1] = -1, 1
   np.all(((P @ G).sum(axis=0) - tfc) < tol)

   # Flux operator Btilde
   Dhat = np.vstack((np.zeros((1,N+1)), D, np.zeros((1,N+1))))
   Btilde = Dhat + (G.T @ P)

   # Laplacian
   L = D @ G

   # Fist mimetic condition
   np.all((L).sum(axis=1) < tol)

   # Second Mimetic Condition (fails)
   tfc = np.zeros(N+2)
   tfc[0], tfc[-1] = -1, 1
   np.all((L.sum(axis=0) - tfc) < tol)


Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
