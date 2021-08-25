import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt

def secante(f, x0, x1, tol, iterMax):
    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    er = []
    err = tol+1
    k = 0
    j = 1
    xk = x0
    xj = x1
    while err> tol and k<iterMax:
        k = k + 1
        j = k-1
        n = sp.N(f1.subs(x,xk))
        u = sp.N(f1.subs(x,xj))
        d = n-u
        a = sp.N(xk - xj)
        xk = xk - a*n/d
        err = abs(f1.subs(x, xk))
        er.append(sp.N(err))
    
    plt.rcParams.update({'font.size':14})
    ejex = np.arange(1,k+1,1)
    fig, graf = plt.subplots()
    graf.plot(ejex,er, 'b', marker='o',markerfacecolor='red',markersize=10)
    plt.show()
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_k)|$')
    graf.set_title('Metodo de Newton-Rapshon (Iteraciones vrs Error)')
    graf.grid(True)
    return [xk,xj, k, err]

f = 'exp(-x^2)-x'
x0 = 0
x1 = 1
tol = 10**-3
iterMax = 1000

y = secante(f, x0, x1, tol, iterMax)
print(y)