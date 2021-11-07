import sys
from sympy import *


# FUNCION SIMPSON COMPUESTO
#   Funcion que se encarga de aproximar una integral definidad por el metodo de simpson compuesto
#     Parametros de entrada
#     funcion: funcion a evaluar, 
#     invervalo: Valores a y b donde se define la integral 
#     puntos: Puntos donde se aproxima la funcion con el metodo
#     Parametros de salida
#     aprox: valor de la aproximacion con el metodo de simpson compuesto
#     error: valor de la cota de error para el metodo

def simpson_compuesto(funcion, intervalo, puntos):

    #Conversion de la funcion de string a simbolica
    x = Symbol('x') #Inicializa "x" como la variable de la funcion a ingresar
    f = sympify(funcion) #Se traduce la funcion tipo string a una aritmetica
    fx = lambdify(x, f, modules=['numpy']) #Se inicializa la funcion

    a = intervalo[0] #Se extrae el valor inicial del intervalo
    b = intervalo[1] #Se extrae el valor final del intervalo
    puntos = 7 #Cantidad de puntos

    h = (b-a)/(puntos-1) #Se calcula el valor de "h"
    lista_x = [] #Se inicializa la lista de los valores de x (x0, x1, x2, ...)
    
    contador = 0 #Inicializa contador para ciclo

    #Se crea la lista con los puntos necesarios para la funcion dada
    while(contador < puntos):
        
        x_i = a + ( contador*h )  #Se calcula x_i
        lista_x+=[x_i] #Se anade x_i a la lista de valores de x para el metodo
        contador+=1 #Aumenta contador

    contador = 1 #Inicializa contador para ciclo
    indices_pares = 0 #Inicializa resultado de sumatoria de indices pares
    indices_impares = 0 #Inicializa resultado de sumatoria de indices impares
            
    #Para el metodo se requieren aquellos termninos x_i pares e impares
    #Este ciclo anade los elementos seun su indice a una lista respectiva 
    #Cada indice se evalua en la funcion original para la integral
    while(contador < puntos -1 ): 

        if( contador %2 == 0): #Verifica que el indice sea par
            indices_pares += fx(lista_x[contador]) #Sumatoria indices pares para el metodo
        else: #Si el indice no es par
            indices_impares += fx(lista_x[contador]) #Sumatoria indices impares para el metodo
        contador += 1 #Aumenta contador

    #Aproximacion de la integral por el metodo de Simpson Compuesto 
    aproximacion = (h/3)*(fx(lista_x[0]) + 2*indices_pares + 4*indices_impares + fx(lista_x[-1])) #Calcula aproximacion para el metodo de Simpson Compuesto
    
      #Seccion de derivadas para obtener la cuarta funcion para el error
    derivada_1 = f.diff(x) #Se calcula la primera derivada de f
    derivada_2 = derivada_1.diff(x) #Se calcula la segunda derivada de f
    derivada_3 = derivada_2.diff(x) #Se calcula la tercera derivada de f
    derivada_4 = derivada_3.diff(x) #Se calcula la cuarta derivada de f
    derivada_4_x = lambdify(x, derivada_4, modules=['numpy']) #Se inicializa la funcion de la cuarta derivada de f
    
    #Calculo del error para el metodo de de Simpson Compuesto
    error = (((b - a) * (h**4))/180) * abs(derivada_4_x(2)) #Calcula error

    print([aproximacion, error]) #Muestra en pantalla la aproximacion y el error
    return[aproximacion, error] #Retorna aproximacion y error

#Funcion de prueba 1
simpson_compuesto('log(x)', [2,5], 7) #Ejemplo

#Funcion de prueba 2
simpson_compuesto('sin(x) / x', [1,2], 11) #Ejemplo