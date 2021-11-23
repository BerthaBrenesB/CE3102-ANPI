from sympy import sympify, Symbol, lambdify, Abs
from numpy import linspace
import numpy as np
import sympy as sp 
from sympy import *
import matplotlib.pyplot as plt

# Metodo de Runge-Kutta para aproximar la solucion de un problema de valor inicial
# a,b ->Valor del rango donde se evaluara el metodo
# y0 -> valor incial de y 
# m => los puntos del intervalo
# La salida es Pares -> Pares x y
def runge_kutta_4(a,b,y0,m):
    x = Symbol('x') #Inicializa "x" como la variable de la funcion a ingresar
    f = sympify('(x+y)/x') #Se traduce la funcion tipo string a una aritmetica
    h = (b-a)/(m-1)
    x_v = linspace(a,b,int(m)+1) # crear un vector column con u valor inicial a y un vector final b cada h espacio
    
    y_v = [y0] # creo el vecto de y con el valor inical
    for n in range(0,int(m)):
            k1 = f.subs(x,x_v[n]) # inicializo el valor de K1 y en evualo en x
            k1 = k1.subs('y', y_v[n]) # Evalueo en y
            k2 = f.subs(x,x_v[n]+ h/2) # inicializo el valor de K2 y en evualo en x
            k2 = k2.subs('y',y_v[n] + h*(k1/2)) # Evalueo en y
            k3 = f.subs(x,x_v[n]+ h/2) # inicializo el valor de K3 y en evualo en x
            k3 = k3.subs('y',y_v[n]+ h*(k2/2)) # Evalueo en y
            k4 = f.subs(x,x_v[n]+ h) # inicializo el valor de K4 y en evualo en x
            k4 = k4.subs('y', y_v[n] + h*k3) # Evalueo en y
            y = y_v[n] +(h/6)*(k1 + 2*k2 + 2*k3 + k4); # determino y evaluo en x el valor de yk+1
            y_v.append(y);# concateno el valor a y
    plt.rcParams.update({'font.size':14})
    plt.plot(x_v,y_v, marker='o',color='red')
    plt.xlabel('Puntos en el intervalo')
    plt.ylabel('Puntos de la imagen que aproximan la solución')
    plt.show()
    return [y_v,x_v]

print(runge_kutta_4(2,10,4,50))
