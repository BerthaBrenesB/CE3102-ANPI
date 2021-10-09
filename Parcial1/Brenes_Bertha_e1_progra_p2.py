import numpy
import numpy as np
from sympy.matrices.expressions.determinant import det


# Calculo de la factorizacion L y U de A
def fact_LU(A):
    n = A.shape[0]#Obtiene el numero de filas de A
    #Obtiene la matriz de inferior y superior en version LU
    U=A
    L=np.zeros((n,n))
    for k in range(0,n-1):#Columnas
        for i in range(k+1,n):#Filas
            M_i_k=U[i,k]/U[k,k]
            L[i,k]=M_i_k#Agrega a la matriz los valores ocupuesto a los "Multiplicadores" de la matriz
            for j in range(k,n):#Ciclo para agregar los nuevos valores de la matriz U
                U[i,j]=U[i,j]-M_i_k*U[k,j]
                if i==j:
                    L[i,j]=1
    L[0,0]=1
    return [L, U]
# Calculo del determinante de la matriz A
def det_fact_lu(A):
    m = A.shape[1]
    det_A = 1 # Variable neutra para la multiplicacion
    for j in range(0,m-1):
        det_A = A[j,j]*det_A; # iteracion de A(1,1)*..*A(m,m)
    return det_A;
    
        




#Ejemplo matriz 5x5
A = np.matrix('25 15 -5 -10 -11;15 10 1 -7 -6;-5 1 21 4 3;-10 -7 4 18 17;-10 -7 4 18 17')
print("Calculo de la factorizacion:")
[L,U] = fact_LU(A)
print(L)
print(U)
print("Calculo del determinante:")
d = det_fact_lu(A)
print(d)
