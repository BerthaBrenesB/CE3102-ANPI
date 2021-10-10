import math
import sympy as sp 
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

#Parametros de entrada
#A matriz nxn, #b matriz de valores independientes
#Salida
#x matriz  de las soluciones
def fact_lu(A, b):
    n = A.shape[0]#Obtiene el numero de filas de A
    m = A.shape[1]#Obtiene el numero de columnas de A
    t = b.shape[0]#Obtiene el numero de filas de b
    o = b.shape[1]#Obtiene el numero de columnas de B
    if n != m:#Comprueba que la matriz sea cuadrada
        raise ValueError("La matriz no es cuadrada")
    if t != n or o!=1:#Comprueba que las matriz b sea las dimensiones adecuadas
        raise ValueError("La matriz b no tiene las dimensiones correctas")
    if teorema_2_lu(A,n) == False:#Comprueba las submatrices sean invertibles
        raise ValueError("La matriz no cumple el teorema 2 de Fact LU")

    M = matriz_inf_lu(A, n)#Obtiene la matriz de inferior y superior en version LU

    y = hacia_Adelante(M[1], b, n)#Resuelve Ly=b
    x = hacia_Atras(M[0], y, n)#Resuleve Ux=y
    return x

#Parametros de entrada
#A matriz  cuadrada, #n largo de la matriz
#Salida
#Valor booleado que identifica si se cumple o no el teorema 2
#Evaluacion del Teorema de LU
def teorema_2_lu(A, n):#Comprueba que las submatrices tiene determinante diferente de 0
    Valor = True#Booleano que sirve para identificar si se cumplio o no el teorema
    for k in range(0, n):
        if np.linalg.det(A[0:k, 0:k]) == 0:
            Valor = False#Se detiene si una submatriz no cumple el teorema 2
            break
    return Valor#True se cumple, #False se incumple


#Parametros de entrada
#A matriz cuadrada, #n largo de la matriz
#Salida
#U matriz Triangular superior LU, #L Matriz Triangular inferior LU
#Matrices inferior y superior version Factorizacion de LU
def matriz_inf_lu(A,n):
    U=A
    L=np.zeros((n,n))

    for k in range(0,n-1):#Columnas
        for i in range(k+1,n):#Filas
            Mik=U[i,k]/U[k,k]
            L[i,k]=Mik#Agrega a la matriz los valores ocupuesto a los "Multiplicadores" de la matriz
            for j in range(k,n):#Ciclo para agregar los nuevos valores de la matriz U
                U[i,j]=U[i,j]-Mik*U[k,j]
                if i==j:
                    L[i,j]=1
    L[0,0]=1
    return U,L


#Parametros de entrada
#A matriz  triangular superior, #b nueva matriz de valores independientes, #n largo de la matriz
#Salida
#soluciones matriz de soluciones "x#"
#Sustitucion hacia atras
def hacia_Atras(A, b, n):
    soluciones = np.zeros((n, 1))
    i = n - 1;

    while (i >= 0):
        sumatoria = 0
        for j in range(i + 1, n):
            sumatoria = sumatoria + (A[i, j] * soluciones[j])

        x = (1 / A[i, i]) * (b[i] - sumatoria)
        soluciones[i] = x
        i = i - 1
    return soluciones


#Parametros de entrada
#A matriz  triangular inferior, #b nueva matriz de valores independientes, #n largo de la matriz
#Salida
#soluciones matriz de soluciones "x#"
#Sustitucion hacia adelante
def hacia_Adelante(A, b, n):
    soluciones = np.zeros((n, 1))#Matriz de soluciones
    soluciones[0] = b[0] / A[0, 0]#Primera solucion
    i = 1;

    while (i < n):#Mueve las filas
        sumatoria = 0
        for j in range(0, i):#Mueve columnas
            sumatoria = sumatoria + (A[i, j] * soluciones[j])#Sumatoria necesaria

        x = (b[i] - sumatoria) / A[i, i]#Calcular el valor de la solucion "x"
        soluciones[i] = x#Agrega la solucion a la matriz de soluciones
        i = i + 1
    return soluciones


# Calculo de xk por medio del metodo de Newton Rapshon
# x0 es un vector con los valores iniciales
# f es un vector de funciones 
# v es un vector of variables
# tol is the maximum tolerance allowed for the error
# iterMax la cantidad de iteraciones maxima
def newton_raphson(x0,f,v,tol,iterMax):
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    z = sp.Symbol('z')
    xk = np.transpose(np.matrix(x0)) # convertir el vector en una matriz fila
    n = len(f) # cantidad de funciones y numero de filas
    m = len(v) # Cantidad de variables y numero de columnas
    df = [] # variable para las funciones diferenciales
    err = [] # variable para los errores
    f2 = []
    if(len(x0) != len(v)):
        raise ValueError('No coincide la cantidad de variables con el vector inicial')
    for i in range(0,n): # iteraciones con respecto a las funciones
        f1 = sp.sympify(f[i]) # Conversion de cada funcion
        f2.append(f1)
        rows = [] #Crear las lista de filas
        for j in range(0,m): # iteracion de variables
            g = sp.diff(f1,v[j]) # Calculo de la derivada de la funcion 
            rows.append(g)
        df.append(rows) # creacion de la matriz con derivadas
    for item in range(0,iterMax): # Calculo de iteraciones
        fx = [] # variable para el valor de fx
        print(item)
        for i in range(0,n):  
            f1 = sp.sympify(f[i])
            for j in range(0,m): # Variables
                f1 = f1.subs(v[j],xk[j]) # Calculo de fx con la evaluacion de x0 con respecto a la variable

            fx.append(float(f1))
        Jk = np.zeros((n,m)) # Matriz de jacobiana
        for j in range(0,len(df)):
            for i in range(0,m):
                #print(df[j][i])
                #print(type(df[j][i]) )
                if(type(df[j][i]) != int): # se calcula si el val9or de la funcion derivada es un entero y por lo tanto no se evalua
                    Jk[j,i] = df[j][i].subs(v[i],xk[i]) # se evalua en la diferencial el valor inicial
                else:
                    Jk[j,i] = df[j][i] # se coloca el valor de la derivada 
        fx = np.matrix(fx)
        b = np.transpose(fx)
        y = fact_lu(Jk,b) # Calculo de y por medio del metodo LU
        xk = xk - y
        print(xk)
        error = np.linalg.norm(np.linalg.norm(fx)) # calculo del error
        err.append(error)
        if(np.linalg.norm(fx)< tol): # Calculo de la tolerancia del error
            break
        
    # Grafica del error contra las iteraciones
    plt.rcParams.update({'font.size':14})
    ejex = np.arange(0, item+1)
    plt.plot(ejex,err, marker='o',color='red')
    plt.xlabel('Iteraciones ($k$)')
    plt.ylabel('$|f(x_k)|$')
    plt.title('Metodo de Newton Raphson Iteraciones vrs Error)')
    plt.show()
    return [xk,item,error]

A = ('3*x**2 -4*y + z**2','x**2 +y**2+z**2 -1','2*x**2 + y**2 -4*z')
x0 = (0.5,0.5,0.5)
x = sp.Symbol('x')
y = sp.Symbol('y')
z = sp.Symbol('z')
v = (x,y,z)
print(newton_raphson(x0,A,v,10**-5,10))


tol = 10**(-10)
itermax = 500
var1 = ["x1","x2"]
f1 = ["exp(x1^2)-exp(sqrt(2)*x1)","x1-x2"]
f1x0 =(2.3,2.3)
#print(newton_raphson(f1x0,f1,var1,tol,itermax))


f2 = ["x1+exp(x2)-cos(x2)","3*x1-x2-sin(x2)"]
f2x0 = (1.5,2)
#print(newton_raphson(f2x0,f2,var1,tol,itermax))


f3 = ["x1^2-2*x1-x2+0.5","x1^2+4*x2^2-4"]
f3x0 = (3,2)
#print(newton_raphson(f3x0,f3,var1,tol,itermax))


f4 = ["x1^2+x2^2-1","x1^2-x2^2+0.5"]
f4x0 = (0.7,1.2)
#print(newton_raphson(f4x0,f4,var1,tol,itermax))


f5 = ["sin(x1)+x2*cos(x1)","x1-x2"]
f5x0 = (1.2,-1.5)
#print(newton_raphson(f5x0,f5,var1,tol,itermax))

f6 = ["x2*x3+x4*(x2+x3)","x1*x3+x4*(x1+x3)","x1*x2+x4*(x1+x2)","x1*x2+x1*x3+x2*x3-1"]
f6x0 = [-1,-1,-1,-1]
var6 = ["x1","x2","x3","x4"]
#print(newton_raphson(f6x0,f6,var6,tol,itermax))

