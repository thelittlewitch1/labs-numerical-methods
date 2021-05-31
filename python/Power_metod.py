import numpy as np
from numpy import linalg as la


def fun(M, u0, e):
    lambda0 = 0
    lambda1 = max(np.dot(M, u0))
    k = 0
    while abs(lambda0 - lambda1) > e:
        lambda0 = lambda1
        k = k + 1
        lambda1 = max(np.dot(la.matrix_power(M, k + 1), u0)) / max(np.dot(la.matrix_power(M, k), u0))
    return np.dot(la.matrix_power(M, k + 1), u0), lambda1


def isSA(N):
    flag = True
    if (N.ndim == 2) and (N.shape[0] == N.shape[1]):
        n = N.shape[0]
        for i in range(0, n):
            for j in range(0, n):
                if N[j, i] != N[i, j].conjugate():
                    flag = False
    else:
        flag = False
    return flag


if __name__ == '__main__':
    A = np.array([
        [4, 1, 1],
        [1, 2, 1],
        [1, 1, 3],
    ])
    u = np.array([[1], [0], [1]])
    eps = 10 ** -3

    if isSA(A) and (A.shape[1] == u.shape[0]):
        print("A:\n", A, "\nu:", u)
        [u1, lambd] = fun(A, u, eps)
        print("total max lambda:", lambd, '\ncorresponding vec: \n', u1)

        print("A^-1:\n", la.inv(A), "\nu:", u)
        [u1, lambd] = fun(la.inv(A), u, eps)
        print("total min lambda:", lambd, '\ncorresponding vec: \n', u1)
    else:
        print('Matrix A: \n', A, '\nis not self-adjoint (Hermitian) matrix, or something else is wrong.')
