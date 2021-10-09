import sympy as sp 
from sympy import *
import math
from math import e

def metodo_biseccion_mod(a,b,n,tol):
    x = sp.Symbol('x')
    f1 = 2 + cos(e**x -2) -e**x
    xk = 0
    h = 0
    if (sp.N(f1.subs(x,a))*sp.N(f1.subs(x,b)) < 0): #Se comprueba la condicion de Bolzano
        print('nuevo intervalo')
        for j in range(0,n-1):
            h = (b-a)/n # se calcula el h
            c_j = a + j*h # Se calcula el c_j
            c_j1 = a + (j+1)*h # Se calcula el c_j+1
            if sp.N(f1.subs(x,c_j))*sp.N(f1.subs(x,c_j1)) < 0: # se comprueba si el cero de f se encuentra si no se cumple se continua con nuevo j
                xk = (c_j1+c_j)/2 # Se calcula el xk
                a = c_j # Se actualiza el valor de a
                b = c_j1 # Se actualiza el valor de b
            error = abs(sp.N(f1.subs(x,xk)))
            if error <tol:
                return xk
        return xk
    else: # El intervalo no cumple la funcion de bolsano
        xk = 'NA'
        print('Intervalo seleccionado no cumple con la funcion de bolsano')
        return xk

x = metodo_biseccion_mod(0,1.2,10,10**-5)
print(x)
