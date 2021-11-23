import numpy
import numpy as np
import sys
from sympy import *
import sympy as sp
from scipy import optimize

#Funcion que calcula la cota de error del trazadir cubico
#Entradas: Una funcion, y el número de puntos
#Salidas: Cota de erro

def cota_traz_cubico(fentrada, xv):
    # Cambiamos la funcion a simbolico
    fs = sp.sympify(fentrada)
    print("Funcion a utilizar: " , end ="")
    print(fs)
    print("Con puntos: ", end ="")
    print(xv)
    #El primer paso es calcular las distancia máxima entre puntos
    #Calucaremos todas las distancias y luego obtenemos el mayor
    n = len(xv)  #Largo del vector de puntos
    dist = []     #Aqui se guardarán las distancias
    #For para calcular las distancias
    for i in range (0, n-1):
        dist.append(xv[i+1] - xv[i])
    h = max(dist) #Valor maximo de distancias = h

    #Calcularemos la cuarta derivada de la funcion de entrada
    fs4 = sp.diff(fs,'x', 4)
    #Buscamos valores extremos del intervalo 
    a = xv[0]
    b = xv[n-1]
    
    faux = -(abs(fs4))  #Creamos una funcion auxiliar negativa de la derivada
    fnumeric = sp.lambdify('x',faux) #Pasamos la funcion a algo que entienda Scipy
    #Calculamos el valor en x donde la funcion derivada en max
    x_max = optimize.fminbound(fnumeric, a, b)

    #Pasamos la funcion derivada a numeral
    fs4aux = sp.lambdify('x',abs(fs4))

    #Calculamos la cota de error con las variable calculadas
    cota_del_error = ((5*(h**4))/384)*fs4aux(x_max)
    print("Cota del error: " , end ="")
    print(cota_del_error)
    
funcion = 'exp(x/2)' #Funcion a utilizar
xv = [1, 1.5, 1.75, 2.15, 2.4, 3] #Puntos con los que se va a trabajar
                                  #Intervalo de [1,3]
cota_traz_cubico(funcion, xv)












