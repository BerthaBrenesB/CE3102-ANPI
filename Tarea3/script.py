import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt


def newton_raphson( x0, tol, iterMax):
    x = sp.Symbol('x')
    f1 =  sp.sympify('2.4719*x**2 -432.49*x  + 16450')
    df1 = sp.diff(f1,x)
    er = []
    err = tol+1
    k = 0
    xk = x0
    while err> tol and k<iterMax:
        k = k+1
        n = sp.N(f1.subs(x,xk))
        d = sp.N(df1.subs(x,xk))
        xk = xk-n/d
        err = abs(f1.subs(x,xk))
        er.append(sp.N(err))

    plt.rcParams.update({'font.size':14})
    ejex = np.arange(1,k+1,1)
    fig, graf = plt.subplots()
    graf.plot(ejex,er, 'b', marker='o',markerfacecolor='red',markersize=10)
    
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_k)|$')
    graf.set_title('Metodo de Newton-Rapshon (Iteraciones vrs Error)')
    graf.grid(True)
    plt.show()
    return [xk, k, err]


x0 = 50
tol = 10**-9
iterMax = 1000

y = newton_raphson(x0,tol, iterMax)
print(y)