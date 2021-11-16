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
    xk = np.matrix(xk) # Creo una matriz a partir de los valores de xk
    yk = np.matrix(yk) # Creo una matriz a partir de los valores de yk
    [m1, m2] = xk.shape # Determino el tamano de la matriz
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
            d = yk[j+1]- yk[j] # Se calcula el dividendo 
            b = xk[j+i] - xk[j]
            nuev2 = np.append(nuev2,d/b)
        iter = iter - 1
        poli_inter = poli_inter + nuev2[1]*v
        yk = nuev2
    poli_inter = expand(poli_inter)
    return poli_inter
# Funcion que realiza el metodo de euler
# f : funcion a evaluar, intervalo: Valor del rango donde se evaluara el metodo
# h : funcion de h ya calculada, y0: valor incial de y 
# La salida es Pares: Pares x y, poli: polinomio de interpolacion es tipo symboli
def euler(f, intervalo, h,y0):
    x = Symbol('x') #Inicializa x como la variable de la funcion a ingresar
    f = sympify(f) #Se traduce la funcion tipo string a una aritmetica
    a = intervalo[0] # determino a
    b = intervalo[1] # determino b
    nump = (b-a)/h # Determino el numero de puntos
    xv = linspace(a,b,int(nump)) # crear un vector column con u valor inicial a y un vector final b cada h espacios
    yv = [y0] # creo el vector de y con el valor inical
    for n in range(0,int(nump)-1): 
            y = yv[n] +h*f.subs(x,xv[n]) # determino el calculo de yk+1 y  evaluo x en la funcion
            y = y.subs('y', yv[n])  # evaluo y en la funciono
            yv.append(y); # concateno el valor a y
    poli_inter = dd_newton(xv,yv); # determino el valor del polinomio
    plt.rcParams.update({'font.size':14})
    plt.plot(xv,yv, marker='o',color='red')
    plt.xlabel('Polinomio de interpolacion')
    plt.ylabel('intervalo')
    plt.show() 
    return [yv,xv],poli_inter
            
print(euler('y-x^2+1',[0,5],0.5,0.5))