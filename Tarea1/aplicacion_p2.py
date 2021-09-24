# def newton_G_m1 Halley's method
import math
import sympy as sp 
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

def newton_raphson(f, x0, tol, iterMax):
    """
    La entradas del metodo son:
    - f: Es la funcion de forma simbolica
    - xo: El valor inicial de xk
    - tol: La tolerancia permitida al error en la funcion
    - iterMax: La cantidad de iteraciones maximas permitidas
    """
    f1 = f
    # Determinacion de la derivada de la funcion
    df1 = sp.diff(f1,x)
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
        # Calculo de xk
        xk = xk-n/d
        # Calculo del error
        err = abs(f1.subs(x,xk))
        # agregar el valor del error a la lista
        er.append(sp.N(err))
    return xk

# Calculo de la distancia en redes de sensores inalambricos basado en el articulo Estimating distances via received signal strength and connectivity in
# wireless sensor networks
def estimating_distance(f, x0, tol, iterMax):
    """
    La entradas del metodo son:
    - fun: Es la funcion en string
    - xo: El valor inicial de dk
    - tol: La tolerancia permitida al error en la funcion
    - iterMax: La cantidad de iteraciones maximas permitidas
    """
    print('Se va a calcular la funcion', f)
    print('Con el methodo One point de tercer orden y la funcion de peso Halley')
    f1 =f
    # definicion de parametros por defectos
    r=10
    a = 4
    odB =4
    v =4
    x1 = 7
    x2 =6
    # primer calculo de dk
    dk = x0
    # Calculo de S
    S = math.pi*r**2
    # Calculo de K
    k = 10*a / math.log(10)
    # Calculo de sigma R
    oR = odB**2/ (10*a)**2
    # Calculo iterativo de dk por medio de newton raphson
    dk = newton_raphson(f1,dk,tol,iterMax)
    # Calculo de g n funcion de d
    gd = (2*S /math.pi)* math.acos(dk/(2*r)) -dk* math.sqrt(r**2 - (dk**2/4))
    # calculo de sigma c
    oc = (gd**2/(2*v*k**2)) * (1/gd + 1/S)
    # Calculo de la estimacion
    fd = math.log10(x1/dk)/(math.log(10)*oR**2) + (dk*(x2-dk))/(oc**2)
    # calculo del error
    err = abs(f1.subs(x,dk))
    print('Calculo de dk') 
    print(dk)
    print('Calculo del error')
    print(err)
    print('Calculo de la distancia d')
    print(fd)
    return [fd, err]


x = sp.Symbol('x')
f = cos(x) -x
x0 = 1
tol = 10**-5
iterMax = 100
y = estimating_distance(f,x0,tol,iterMax)
