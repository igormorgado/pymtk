# Gradient order 2
# N = 6
P = np.diag((3/8, 9/8, 1, 1, 1, 9/8, 3/8))
# N+1 x N+2
hG = np.array([
    [-8/3, 3, -1/3,  0,  0,   0,  0,    0],
    [0,   -1,    1,  0,  0,   0,  0,    0],
    [0,    0,   -1,  1,  0,   0,  0,    0],
    [0,    0,    0, -1,  1,   0,  0,    0],
    [0,    0,    0,  0, -1,   1,  0,    0],
    [0,    0,    0,  0,  0,  -1,  1,    0],
    [0,    0,    0,  0,  0, 1/3, -3,  8/3]
    ])

# N x N+1
Q = np.diag(np.ones(6))
hD = np.array([
    [-1,  1,  0,  0,  0,  0,  0],
    [ 0, -1,  1,  0,  0,  0,  0],
    [ 0,  0, -1,  1,  0,  0,  0],
    [ 0,  0,  0, -1,  1,  0,  0],
    [ 0,  0,  0,  0, -1,  1,  0],
    [ 0,  0,  0,  0,  0, -1,  1]
])

hDhat = np.vstack((np.zeros(7), hD, np.zeros(7)))


(P @ hG).sum(axis=0)
(P @ hG).sum(axis=1)
(Q @ hD).sum(axis=0)
(Q @ hD).sum(axis=1)
Btilde = hDhat + (hG.T @ P)
Btilde

