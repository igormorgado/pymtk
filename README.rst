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

How to install
--------------

Just type:

.. code-block:: bash

   pip3 install  "git+https://github.com/igormorgado/pymtk/#egg=pymtk"


And it will install all requirements


How to use
----------

To checkout call `ipython3` and run the code below:

.. code-block:: python

   import numpy as np
   import pymtk


   # Keep things nicer
   np.set_printoptions(suppress=True)
   np.set_printoptions(precision=4)
   np.set_printoptions(linewidth=180)
   
   # Operator Order
   k = 4

   # Discretized points
   N = 4*k-2

   # Anti diagonal matrix
   P_ = lambda N: np.fliplr(np.diag(np.ones(N)))

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
   
   # Is Divergent Skew-Centro-Simmetric ( P_N @ D @ P_{N+1} = - D )
   np.all( (P_(N) @ D @ P_(N+1) ) + D < tol)

   # Gradient

   # Fist mimetic condition
   G_DCZ = (P @ G).sum(axis=1) 
   np.all(G_DCZ < tol)

   # Second Mimetic Condition
   tfc = np.zeros(N+2)
   tfc[0], tfc[-1] = -1, 1
   G_TFC = (P @ G).sum(axis=0) 
   np.all((G_TFC - tfc) < tol)

   # Is Gradient Skew-Centro-Simmetric ( P_{N+1} @ G @ P_{N+2} = -G )
   np.all( (P_(N+1) @ G @ P_(N+2) ) + G < tol)


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

* Once operator is created, for example

.. code-block:: python

   import pymtk
   import numpy as np
   D_4 = pymtk.Divergent(order=4)


Is possible to extract useful operator informations as

1. Upper left(and bottom right) boundary rows

.. code-block:: python

   D_4.boundary_rows
   # OUT
   # array([[-0.915061633,  0.700308166,  0.391050334, -0.224383667,  0.049691834, -0.001605033],
   #        [ 0.041666667, -1.125      ,  1.125      , -0.041666667,  0.         ,  0.         ],
   #        [ 0.         ,  0.041666667, -1.125      ,  1.125      , -0.041666667,  0.         ],
   #        [ 0.         ,  0.         ,  0.041666667, -1.125      ,  1.125      , -0.041666667]])

   - np.flipud(np.fliplr(D_4.boundary_rows))
   # OUT
   # array([[ 0.041666667, -1.125      ,  1.125      , -0.041666667, -0.         , -0.         ],
   #        [-0.         ,  0.041666667, -1.125      ,  1.125      , -0.041666667, -0.         ],
   #        [-0.         , -0.         ,  0.041666667, -1.125      ,  1.125      , -0.041666667],
   #        [ 0.001605033, -0.049691834,  0.224383667, -0.391050334, -0.700308166,  0.915061633]])



2. Inner product weights and associated vector/matrix

.. code-block:: python

   D_4.lambda_
   # OUT
   # array([-0.001808449])
   
   D_4.weights
   # OUT
   # array([1.126736111, 0.744791667, 1.171875   , 0.956597222])

   D_4.weight_vector(11)
   # OUT
   # array([1.126736111, 0.744791667, 1.171875   , 0.956597222, 1.         , 1.,
   #        1.         , 0.956597222, 1.171875   , 0.744791667, 1.126736111])

   np.set_printoptions(precision=5)
   np.diag(D.weight_vector(11))
   # OUT
   # array([[1.12674, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ],
   #        [0.     , 0.74479, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ],
   #        [0.     , 0.     , 1.17187, 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ],
   #        [0.     , 0.     , 0.     , 0.9566 , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ],
   #        [0.     , 0.     , 0.     , 0.     , 1.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     ],
   #        [0.     , 0.     , 0.     , 0.     , 0.     , 1.     , 0.     , 0.     , 0.     , 0.     , 0.     ],
   #        [0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 1.     , 0.     , 0.     , 0.     , 0.     ],
   #        [0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.9566 , 0.     , 0.     , 0.     ],
   #        [0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 1.17187, 0.     , 0.     ],
   #        [0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.74479, 0.     ],
   #        [0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 1.12674]])


3. Operator Vandermonde generators and stencil

.. code-block:: python

   D_4.boundary_generator()
   # OUT
   # array([[-0.5,  0.5,  1.5,  2.5,  3.5,  4.5],
   #        [-1.5, -0.5,  0.5,  1.5,  2.5,  3.5],
   #        [-2.5, -1.5, -0.5,  0.5,  1.5,  2.5],
   #        [-3.5, -2.5, -1.5, -0.5,  0.5,  1.5]])
   
   D_4.stencil
   # OUT
   # array([ 0.041666667, -1.125      ,  1.125      , -0.041666667])


4. Operator Nullspace

.. code-block:: python

   D_4.Nu
   # OUT
   # array([[ -1.,   5., -10.,  10.,  -5.,   1.]])


5. The operator discretized in N intervals

.. code-block:: python

   np.set_printoptions(precision=4)
   D_4(11)
   # OUT
   # array([[-0.9151,  0.7003,  0.3911, -0.2244,  0.0497, -0.0016,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417, -0.    , -0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , -0.    ,  0.0417, -1.125 ,  1.125 , -0.0417, -0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , -0.    , -0.    ,  0.0417, -1.125 ,  1.125 , -0.0417],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0016, -0.0497,  0.2244, -0.3911, -0.7003,  0.9151]])


6. Same for the Gradient (here only the full matrix and weights for reference) 

.. code-block:: python

   import pymtk
   import numpy as np

   np.set_printoptions(precision=9)
   
   G_4 = pymtk.Gradient(order=4)

   G_4.boundary_rows
   # Out
   # array([[-3.361740962,  4.398034398, -1.489045864,  0.552641278, -0.102448227,  0.002559378],
   #        [ 0.152380952, -1.291666667,  1.208333333, -0.075      ,  0.005952381,  0.         ],
   #        [ 0.         ,  0.041666667, -1.125      ,  1.125      , -0.041666667,  0.         ],
   #        [ 0.         ,  0.         ,  0.041666667, -1.125      ,  1.125      , -0.041666667]])

   G_4.weight_vector(11)
   # Out
   # array([0.353298611, 1.231770833, 0.893229167, 1.021701389, 1.         , 1.         , 1.,
   #        1.         , 1.021701389, 0.893229167, 1.231770833, 0.353298611])

   np.set_printoptions(precision=4)
   G_4(11)
   # Out 
   # array([[-3.3617,  4.398 , -1.489 ,  0.5526, -0.1024,  0.0026,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.1524, -1.2917,  1.2083, -0.075 ,  0.006 ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417,  0.    ,  0.    ,  0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.0417, -1.125 ,  1.125 , -0.0417, -0.    , -0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , -0.    ,  0.0417, -1.125 ,  1.125 , -0.0417, -0.    ],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , -0.    , -0.006 ,  0.075 , -1.2083,  1.2917, -0.1524],
   #        [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , -0.0026,  0.1024, -0.5526,  1.489 , -4.398 ,  3.3617]])



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
