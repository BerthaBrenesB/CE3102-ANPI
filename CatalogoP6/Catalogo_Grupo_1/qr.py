from scipy import linalg as la
from numpy import matrix, zeros
import numpy as np
import matplotlib.pyplot as plt
"""
    Metodo QR para encontrar todos los autovalores de una matriz
    :param A: Matriz A
    :return: Matriz con todos los autovalores de la matriz A
    """
def metodo_qr(A):
    # Inicializo los valores iniciales
    tol = 10 ** -6
    A = matrix(A)
    A_ant = matrix(zeros(A.shape))
    [n,m] = A.shape
    iter_max = 100
    U = np.eye(n)
    ks = []
    errors = []
    for k in range(0,iter_max):
        # Calculo de Q y R 
        Q, R = la.qr(A)
        # Almaceno el valor anterior de a
        A_ant = A
        # Calculo el nuevo valor de A
        A = R * Q
        # Calculo el valor de U
        U = U*Q
        # Se verifica la condicion de parada
        norma = la.norm(A_ant - A)
        # agrego a la lista los valores de error
        errors.append(norma)
        ks.append(k)
        # Calculo el error
        if norma < tol:
            break
    plt.rcParams.update({'font.size':14})
    plt.plot(ks,errors, marker='o',color='red')
    plt.ylabel('Iteraciones')
    plt.xlabel('Valores de error')
    plt.show() 
    return A


resultado = metodo_qr([[0, 11, -5], [-2, 17, -7], [-4, 26, -10]])
print(resultado) 