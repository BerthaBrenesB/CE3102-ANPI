import numpy as np
import sys
from sympy import *
# FUNCION SIMPSON
# El objetivo de esta funcion es poder aproximar el valor de una integral
# definida en un intervalo a y b dado    
# Parametros de entrada
# funcion funcion integrable
# intervalo interbalo del funcion    
# Parametro de salida
# fresultado aproximacion de la integral
# error: error de la aproximacion
def simpson(funcion, intervalo):
    x = Symbol('x') 
    f = sympify(funcion) 
    fx =  lambdify(x, f, modules=['numpy'])

    a = intervalo[0] 
    b = intervalo[1] 
    
    x0 = a 
    x1 = (a+b)/2 
    x2 = b 
    
    # Valor h para la formula de Simpson
    h = (b-a)/2 
    
    # Implementacion de la formula de Simpson 
    f_resultado = (h/3)*( fx(x0) + 4*fx(x1) + fx(x2) ) 
    
    # Se procede al calculo del error para la funcion
    df = f.diff(x,4) # Se calcula la cuarta derivada de la funcion inicial
    dfx = lambdify(x, df, modules=['numpy']) # Se iniacilaiza la funcion fx

    f1 = abs(dfx(intervalo[0])) # Valores de los puntos
    f2 = abs(dfx(intervalo[1]))
    # Se procede a garantizar la continuidad de la funcion en todo momento

    if (f1 > f2): #Punto maximo
         max_relativo = [intervalo[0], f1] # Asignacion del punto maximo punto maximo
    else:
         max_relativo = [intervalo[1], f2] # Asignacion del punto maximo punto maximo
    # Se debe validar que la funcion no se indefina en ningun momento
    try:
        punto = []
        solucion = np.solve(f.diff(x,1))
        for i in solucion: # Calcula cual de los resultados es el maximo
            if (abs(dfx(i)) > punto):  # Determina si el absoluto de la deriva es mayor al punto analizado
                punto = [i, abs(dfx(i))] # Se obtiene el punto analizado con el maximo 
        
        if (punto[1] > max_relativo[1]):  # Compara el maximo obtenido en la funcion con el maximo de los extremos
            error = (h**5/90)*abs(dfx(punto[0])) # Se aplica la formula del error para el metodo de Simpson
        else:
            error = (h**5/90)*abs(dfx(max_relativo[0]))   # Se aplica la formula del erro para el metodo de Simpson
        
    except:
        error = (h**5/90)*abs(dfx(max_relativo[0]))  # Se aplica el error para el punto si no es igual al maximo
    
    print(f_resultado, error)
    return(f_resultado, error) # Resultado de la aproximacion a la integral

#Funcion de prueba 1 
simpson('ln(x)', [2,5])
    
#Funcion de prueba 2
simpson('13 / (5*x + 4)', [1 , 2])