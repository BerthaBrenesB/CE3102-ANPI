import numpy as np
import sympy as sp

def muller(f1,x0,x1,x2,tol):
    error= 1e3
    x3 = 0
    while error>tol:
        c = f(x2)
        b = ((x0-x2)**2*(f1(x1) - f1(x2)) - (x1-x2)**2*(f1(x0) - f1(x2))) / ((x0-x1)*(x0-x2)*(x1-x2)) # Calcula y almacena el valor de b
        a = ((x1-x2)*(f1(x0) - f1(x2)) - (x0-x2)*(f1(x1) - f1(x2))) / ((x0-x1)*(x0-x2)*(x1-x2)) # Calcula y almacena el valor de a

        x3 = x2- (2*c)/ (b+np.sign(b)*np.sqrt(b**2 -4*a*c))
        error = abs(x3-x2)
        x0 = x1
        x1 = x2
        x2 = x3
    return x3


f = lambda x: np.sin(x) -x/2
y = muller(f,2,2.2,1.9,1e-5)
print(y)