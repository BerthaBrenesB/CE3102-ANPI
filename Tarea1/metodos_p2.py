import math
import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt

# Implementacion del metodo de newton raphson con B diferente de 0
# especificamente con el metodo de Kanwar-Tomar
def newton_H_m1(f, x0, tol, iterMax):
    """
    La entradas del metodo son:
    - fun: Es la funcion en string
    - xo: El valor inicial de xk
    - tol: La tolerancia permitida al error en la funcion
    - iterMax: La cantidad de iteraciones maximas permitidas
    """
    print('Se va a calcular la funcion', f)
    print('Con el methodo One point y la funcion de peso Kanwar-Tomar')
    x = sp.Symbol('x')
    # Calculo de f simbolico
    f1 = sp.sympify(f)
    # Determinacion de la derivada de la funcion
    df1 = sp.diff(f1,x)
    b = 0.5 # definicion de beta
    # creo una lista para almacenar el error
    er = []
    err = tol+1
    # Variable para las iteraciones
    k = 0
    # Calculo del primer valor de xk
    xk = x0
    # Comienzo de las iteraciones
    while err> tol and k<iterMax:
        k = k+1
        # calculo del nominador
        n = sp.N(f1.subs(x,xk))
        # Calculo del denominador
        d = sp.N(df1.subs(x,xk))
        # Calculo de f1/f'1
        a = n/d
        # calculo de xk
        xk = xk - (1/(1 + b*a))*a
        # Calculo del error con respecto a xk
        err = abs(f1.subs(x,xk))
        er.append(sp.N(err))
    print('Calculo de Xk') 
    print(xk)
    print('Calculo del error')
    print(err)
    plt.rcParams.update({'font.size':14})
    ejex = np.arange(1,k+1,1)
    fig, graf = plt.subplots()
    graf.plot(ejex,er, 'b', marker='o',markerfacecolor='red',markersize=10)
    plt.show()
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_k)|$')
    graf.set_title('Metodo de Newton-Rapshon (Iteraciones vrs Error)')
    graf.grid(True)
    return [xk, k, err]


f = 'cos(x) -x'
x0 = 2.1
tol = 10**-5
iterMax = 500


# Implementacion del metodo de newton raphson con B diferente de 0
# especificamente con el metodo de Noor
def newton_H_m2(f, x0, tol, iterMax):
    """
    La entradas del metodo son:
    - fun: Es la funcion en string
    - xo: El valor inicial de xk
    - tol: La tolerancia permitida al error en la funcion
    - iterMax: La cantidad de iteraciones maximas permitidas
    """
    print('Se va a calcular la funcion', f)
    print('Con el methodo One point y la funcion de peso Noor')
    x = sp.Symbol('x')
    # Calculo de f simbolico
    f1 = sp.sympify(f)
    # Determino la derivada de la funcion
    df1 = sp.diff(f1,x)
    # creo una lista para almacenar el error
    er = []
    err = tol+1
    # Variable para las iteraciones
    k = 0
    b = 0.5 # definicion de beta
    xk = x0
    # Comienzo de las iteraciones
    while err> tol and k<iterMax:
        k = k+1
        # calculo del nominador
        n = sp.N(f1.subs(x,xk))
        # Calculo del denominador
        d = sp.N(df1.subs(x,xk))
        u = n/d
        # Calculo de xk
        xk = xk - u*1
        # Calculo del error
        err = abs(f1.subs(x,xk))
        er.append(sp.N(err))

    print('Calculo de Xk') 
    print(xk)
    print('Calculo del error')
    print(err)
    plt.rcParams.update({'font.size':14})
    ejex = np.arange(1,k+1,1)
    fig, graf = plt.subplots()
    graf.plot(ejex,er, 'b', marker='o',markerfacecolor='red',markersize=10)
    plt.show()
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_k)|$')
    graf.set_title('Metodo de Newton-Rapshon (Iteraciones vrs Error)')
    graf.grid(True)
    return [xk, k, err]




# def newton_G_m1 Halley's method
def newton_G_m1(f, x0, tol, iterMax):
    """
    La entradas del metodo son:
    - fun: Es la funcion en string
    - xo: El valor inicial de xk
    - tol: La tolerancia permitida al error en la funcion
    - iterMax: La cantidad de iteraciones maximas permitidas
    """
    print('Se va a calcular la funcion', f)
    print('Con el methodo One point de tercer orden y la funcion de peso Halley')
    x = sp.Symbol('x')
    # Calculo de f simbolico
    f1 = sp.sympify(f)
    # Determino la primera derivada de la funcion
    df1 = sp.diff(f1,x)
    # Determino la segunda derivada de la funcion
    d2f1 = sp.diff(f1,x,x)
    # creo una lista para almacenar el error
    er = []
    err = tol+1
    # Variable para las iteraciones
    k = 0
    # Calculo del primer valor de xk
    xk = x0
    # Comienzo de las iteraciones
    while err> tol and k<iterMax:
        k = k+1
        # calculo del nominador
        n = sp.N(f1.subs(x,xk))
        # Calculo de la primera derivada con respecto a xk
        d = sp.N(df1.subs(x,xk))
        # Calculo de la segunda derivada con respecto a xk
        d2 = sp.N(d2f1.subs(x,xk))
        u = n/d
        # Calculo del peso w
        w = (n*d2)/(d**2)
        # Calculo de xk
        xk = xk - 2*u/(2-w)
        # Calculo del error
        err = abs(f1.subs(x,xk))
        er.append(sp.N(err))

    print('Calculo de Xk') 
    print(xk)
    print('Calculo del error')
    print(err)
    plt.rcParams.update({'font.size':14})
    ejex = np.arange(1,k+1,1)
    fig, graf = plt.subplots()
    graf.plot(ejex,er, 'b', marker='o',markerfacecolor='red',markersize=10)
    plt.show()
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_k)|$')
    graf.set_title('Metodo de Newton-Rapshon (Iteraciones vrs Error)')
    graf.grid(True)
    return [xk, k, err]


# def newton_G_m2 Chebyshev's
def newton_G_m2(f, x0, tol, iterMax):
    """
    La entradas del metodo son:
    - fun: Es la funcion en string
    - xo: El valor inicial de xk
    - tol: La tolerancia permitida al error en la funcion
    - iterMax: La cantidad de iteraciones maximas permitidas
    """
    print('Se va a calcular la funcion', f)
    print('Con el methodo One point de tercer orden y la funcion de peso Chebyshev')
    x = sp.Symbol('x')
    # Calculo de f simbolico
    f1 = sp.sympify(f)
    # Determino la primera derivada de la funcion
    df1 = sp.diff(f1,x)
    # Determino la segunda derivada de la funcion
    d2f1 = sp.diff(f1,x,x)
    # creo una lista para almacenar el error
    er = []
    err = tol+1
    # Variable para las iteraciones
    k = 0
    # Calculo del primer valor de xk
    xk = x0
    # Comienzo de las iteraciones
    while err> tol and k<iterMax:
        k = k+1
        # calculo del nominador
        n = sp.N(f1.subs(x,xk))
        # Calculo de la primera derivada con respecto a xk
        d = sp.N(df1.subs(x,xk))
        # Calculo de la segunda derivada con respecto a xk
        d2 = sp.N(d2f1.subs(x,xk))
        u = n/d
        # Calculo del peso w
        w = (n*d2)/(d**2)
        # Calculo de xk
        xk = xk - u*(1 + (1/2)*w)
        # Calculo del error
        err = abs(f1.subs(x,xk))
        er.append(sp.N(err))

    print('Calculo de Xk') 
    print(xk)
    print('Calculo del error')
    print(err)
    plt.rcParams.update({'font.size':14})
    ejex = np.arange(1,k+1,1)
    fig, graf = plt.subplots()
    graf.plot(ejex,er, 'b', marker='o',markerfacecolor='red',markersize=10)
    plt.show()
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x_k)|$')
    graf.set_title('Metodo de Newton-Rapshon (Iteraciones vrs Error)')
    graf.grid(True)
    return [xk, k, err]

y = newton_G_m2(f,x0,tol, iterMax)
y1 = newton_G_m1(f,x0,tol, iterMax)
y2 = newton_H_m1(f,x0,tol, iterMax)
y2 = newton_H_m2(f,x0,tol, iterMax)