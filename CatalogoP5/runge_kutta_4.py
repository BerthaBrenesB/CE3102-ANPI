from sympy import sympify, Symbol, lambdify, Abs
from numpy import linspace
import numpy as np
import sympy as sp 
from sympy import *
import matplotlib.pyplot as plt
# Funcion que realiza el metodo de Diferencias Divididas de Newton
# puntos:una matriz mx2 donde la column 1 son los valores de x 
# Y la columna 2 son los valores de y
def dd_newton(xk, yk):
    xk = np.matrix(xk)
    yk = np.matrix(yk)
    [m1, m2] = xk.shape
    [n1,n2] = yk.shape
    if m1 != n1: # Comprueba las lista de pares ordenados
        raise ValueError('No son pares ordenados')
    x = sp.Symbol('x')
    puntos = np.append(xk,yk,axis=1)
    poli_inter = puntos[0,1] # Se almacena la primera diferencia dividida
    v = 1
    iter = n1-1
    for i in range(1,n1):
        v = v*(x- puntos[i-1,0]) # Se calcula la variable que sera multiplicado por las diferencias
        nuev2 = np.zeros(1) # Se almacenara los multiplicadores d/b
        for j in range(0,iter):
            d = yk[j+1]- yk[j] # Se calcula el dividendo ejemplo: F[x1,x2]-[F[x0,x1]
            b = xk[j+i] - xk[j]
            nuev2 = np.append(nuev2,d/b)
        iter = iter - 1
        poli_inter = poli_inter + nuev2[1]*v
        yk = nuev2
    poli_inter = expand(poli_inter)
    return poli_inter


def runge_kutta_4(f,intervalo, h):
    x = Symbol('x') #Inicializa "x" como la variable de la funcion a ingresar
    f = sympify(f) #Se traduce la funcion tipo string a una aritmetica
    a = intervalo[0]
    b = intervalo[1]
    nump = (b-a)/h
    xv = linspace(a,b,int(nump)+1) # crear un vector column con u valor inicial a y un vector final b cada h espacios
    yv = [0.5]
    for n in range(0,int(nump)):
            k1 = f.subs(x,xv[n])
            k1 = k1.subs('y', yv[n])
            k2 = f.subs(x,xv[n]+ h/2)
            k2 = k2.subs('y',yv[n] + h*(k1/2))
            k3 = f.subs(x,xv[n]+ h/2)
            k3 = k3.subs('y',yv[n]+ h*(k2/2))
            k4 = f.subs(x,xv[n]+ h)
            k4 = k4.subs('y', yv[n] + h*k3)
            y = yv[n] +(h/6)*(k1 + 2*k2 + 2*k3 + k4);
            yv.append(y);
    
    poli_inter = dd_newton(xv,yv);
    plt.rcParams.update({'font.size':14})
    plt.plot(xv,yv, marker='o',color='red')
    plt.xlabel('Polinomio de interpolacion')
    plt.ylabel('intervalo')
    plt.show()
    return [yv,xv], poli_inter

print(runge_kutta_4('-x*y + 4*x/y',[0,1],0.1))