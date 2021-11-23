from scipy import linalg as la
from numpy import matrix
import numpy as np
from sympy.solvers.diophantine import norm
### Metodo de rayleigh que calcula el valor propio mas pequeño de una matriz
def rayleigh(A):
    # valor de tolerancia
    tol = 10 ** -6
    # conversion a una matriz A
    A = matrix(A)
    [n,m] = A.shape
    # creacion de un vector x
    xk = np.ones(n)
    iter_max = 100
    Id = np.eye(n)
    ok = 0;
    norma = 1
    while norma>tol:
        # Calculo de la transpuesta
        xk_t = np.transpose(xk)
        d1 = np.dot(xk_t,A)
        d = np.dot(d1,xk);
        n = np.dot(xk_t,xk)
        # calculo sigma ma anterior
        ok_ant = ok;
        # calculo de sigma
        ok = d/n
        ok = ok.item((0,0))
        # Calculo de A - sigma*matrizIdentidad
        V = A - ok*Id;
        # Calculo de yk
        yk = np.linalg.solve(V,xk);
        # Calculo de xk
        xk = yk/ np.linalg.norm(yk,2)
        # calculo de la norma
        norma = la.norm(ok_ant - ok)
        print(norma)
        #Calculo el error
        if norma < tol:
            break
    return xk;

A = [[8,1,0,0,0], [1,8,1,0,0], [0,1,8,1,0], [0,0,1,8,1],[0,0,0,1,8]]
print(rayleigh(A))