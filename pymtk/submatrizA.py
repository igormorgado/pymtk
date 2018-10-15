# Algoritmo 9 - Calcula a submatriz A(alpha)


def MatrizA(PiTilde, Nu, q, lambda_):

    k = len(q)
    A = PiTilde.T[:k]

    alpha = lambda_ / q[:len(lambda_)]
    for i in range(len(alpha)):
        A[i] += alpha[i] * Nu[i]

    return A
