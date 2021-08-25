# Implementación de biblioteca, indispensables para correr el programa 
# Por lo tanton
# 1. Se tiene que tener instalado python 3 en la compu
# Tener instalado los paquetes de numpy y matplotlib
# correr el archivo como python3 bisección 
import sympy as sp 
import numpy as np
import matplotlib.pyplot as plt

# Se definen las entradas que en este caso
# tol = es la entrada correspondiente a la tolerancia que tendra el sistema
# f = es la funcion definida en un string
# x0 = es el valor inicial en 0
# iterMax = la cantidad de iteraciones máxima que el sistema realizara 
# La salida esperada de este sistema 
# xk = El punto donde f(xk) = 0 
# e = El valor del error
# Ademas de la funcion ploteada de la cantidad de iteración y el valor de error
def biseccion(f, a, b, tol, iterMax):
# Esta funcion nos va a dar la aproximación de soluciones a una ecuacion 
# f(x) usando el método de la bisección
    x = sp.Symbol('x') # Define del texto x a simbólico
    f1 = sp.sympify(f) # Método Symbol que convierte el texto a simbolico
    er = [] # Array de valores de error
    xk = 0 # Variable que va a contener el calculo de valor de x
    k = 0 # Variable que va a las iteraciones del método
    print('primera valor de f(a)',sp.N(f1.subs(x,a)))
    print('primera valor de f(b)',sp.N(f1.subs(x,b)))
    if sp.N(f1.subs(x,a))*sp.N(f1.subs(x,b)) < 0: # Constantacion de la funcion de bolsano 
        while k<iterMax: # Comienzo del ciclo de iteraciones con parada mientras k sea menor al valor de iteraciones
            k = k+1 # contador de k
            xk=(a+b)/2 # Calculo de xk 
            err = (b-a) / 2**k # Calculo del error
            if sp.N(f1.subs(x,a))*sp.N(f1.subs(x,xk))<0: # Calculo de la mitad mas cercana al cero
                b=xk 
            else:
                a = xk
            er.append(sp.N(err)) # agregar a la lista de errores
            if abs(err)<tol: # Punto de parada cuando el valor del error supera la tolerancia
                break
        plt.rcParams.update({'font.size':14})
        ejex = np.arange(1, k+1)
        plt.plot(ejex,er, marker='o',color='red')
        plt.xlabel('Iteraciones ($k$)')
        plt.ylabel('$|f(x_k)|$')
        plt.title('Metodo de Biseccio Iteraciones vrs Error)')
        plt.show()
        return [xk,err]
    else: # Al no cumplir la funcion de bolsano
        xk = 'NA'
        a ='NA'
        b = 'NA'
        err ="NA"
        print('Intervalo seleccionado no cumple con la funcion de bolsano')
# Calculo de la funcion f(x)=exp(x)-x-2
f = 'exp(x)-x-2'
a = 0
b = 2
tol = 10^-4
iterMax = 100
y = biseccion(f, a,b,tol,iterMax)
print(y)