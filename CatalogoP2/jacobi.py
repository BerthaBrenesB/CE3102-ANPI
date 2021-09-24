import numpy as np

def domin_matri(A):
    n = len(A)
    for f in range(0,n):
        d = A[f,f]
        f1 = A[f,:]
        if abs(d) <= abs(f1).sum():
            return True
        else: 
            return False
        
def jacobi(A,b,xo,tol,iterMax):
    [n,m] = A.shape
    if(n != m):
        raise ValueError('La matriz no es cuadrada')
    t = len(b)
    o = len(b[0])
    if(t !=n and o!=1):
        raise ValueError('La matriz de valores independientes no cumple con las dimensiones de la matriz A')
    if(domin_matri(A) != True):
        raise ValueError('La matriz no es dominante')
    d = np.diag(A)
    D_inv = np.diag(1./d)
    LpU = A - np.diag(d)
    z = D_inv * b.transpose()
    xk = xo
    error = 0
    for k in range(1,iterMax):
        xk = (-D_inv*LpU*xk) +z
        error = np.linalg.norm(A*xk -b)
        if(error<tol):
            break

    return [xk,error]


A = np.matrix('10 1 -1 5;0 -7 1 1; 1 2 5 0;4 4 1 10')
b = np.matrix('1 1 1 1')
iterMax = 150
tol = 1^-15
x0 = np.zeros((4,1))
print(jacobi(A,b,x0,tol,iterMax))