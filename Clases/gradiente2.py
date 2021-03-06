from numpy import linalg, array, absolute
from sympy import sympify, Symbol, diff
from matplotlib import pyplot as plt

# Funcion principal
def gradiente(func, var, xK, tol, iterMax):

    #Función de prueba: (x-2)**4 + (x-2*y)**2", ['x', 'y'], [0, 3], 10**(-10), 13

    if len(var) != len(xK):
        return "La lista de variables y de valores iniciales deben ser iguales."
    
    if tol < 0:
        return "La tolerancia dada debe ser mayor a cero"
    
    funtion = sympify(func) #Se traduce el string "func" a una función de sympy
    itr = 1

    D = []
    A = []

    # Se crea la lista de valores simbolicos
    symbolicList = []
    for i in var:
        symbolicList.append(Symbol(i))
    
    # Se calcula el gradiente de la funcion
    gradient = []
    for i in symbolicList:
        gradient.append(diff(funtion, i))
    
    # Calculo de los valores iniciales
    gk = valueGradient(gradient, var, xK)
    print(gk)
    dk = []
    for i in gk:
        dk.append(i*-1)

    # Se calcula la norma para el primer valor
    norma = linalg.norm(array(gk, dtype='float'), 2)
    
    D.append(itr-1)
    A.append(round(norma, 4))

    while itr <= iterMax:
        # Se calcula el alpha
        ak = calc_ak(funtion, var, xK, dk, gk)

        # Se calcula Xk+1 = xk + a * dk
        ad_k = []
        for i in dk:
            ad_k.append(i * ak)
        xk_1 = []
        for (i, j) in zip(xK, ad_k):
            xk_1.append(i + j)

        # Se calcula el nuevo valor del vector gk
        gk_1 = valueGradient(gradient, var, xk_1)

        # Se calcula el vector para la condicion de parada
        stop_xk_1 = valueGradient(gradient, var, xk_1)

        # Se calcula la norma 2 del vector para la condicion de parada
        norma = linalg.norm(array(stop_xk_1, dtype='float'), 2)

        D.append(itr)
        A.append(round(norma, 4))

        # Se verifica la condicion de parada
        if norma < tol:
            break

        # Se calcula el valor de beta
        bk = calc_bk(gk_1, gk, dk)

        # Se calcula el nuevo valor del vector dk
        bk_dk = [i * bk for i in dk]
        m_gk_1 = [i * -1 for i in gk_1]
        dk = [x1 + x2 for (x1, x2) in zip(m_gk_1, bk_dk)]

        xK = xk_1.copy()
        gk = gk_1.copy()

        itr += 1

    plt.plot(array(D), array(A))
    plt.title("Gráfico")
    plt.legend(["Error"])
    plt.show() 
    
    print(xk_1, round(norma, 4), (itr-1))
    return [xk_1, round(norma, 4)] #Retorna el Xk(aproximación), y el error

# Evalua el gradiente
def valueGradient(gradient, variable, values):

    n = len(variable)
    value = []

    # Se recorre cada una de las derivadas parciales en el gradiente
    for i in range(0, n):
        
        funtion = gradient[i] 

        for j in range(0, n):

            funtion  =  funtion.subs(variable[j], values[j])
        
        value.append(funtion.doit())

    return value    

# Evalua la funcio
def value_funtion(func, var, xK):
    
    n = len(var)
    
    for i in range(0, n):
        func = func.subs(var[i], xK[i])
        print(func)
    return func

# Calcula ak
def calc_ak (func, var, xK, dK, gK):

    # Desigualdad f(xK + aK * dK) - f(xK) <= 0.5 * aK * gK * dK
    print('calc ak')
    print(func, var)
    print(xK, gK, dK)
    ak = 1
    while 1:

        #Calcular primera parte de la desigualdad

        # ak * dk
        ad_K =[]
        for i in dK:
            ad_K.append(i * ak)
        #print(ad_K)
        # xK + aK * dK

        a_ad_K = []
        for (i, j) in zip(xK, ad_K):
            a_ad_K.append(i + j)
        #print(a_ad_K)
        #Evaluar f(xK + aK * dK)
        f_x_ad_k = value_funtion(func, var, a_ad_K)
        #print(f_x_ad_k)
        #Evalar f(xK)
        f_xk = value_funtion(func, var, xK)
        #print(f_xk)
        #Primer parte de la desigualdad
        firstPart = f_x_ad_k - f_xk

        #Calcular segunda parte de la desigualdad

        #Calculo gk * dk
        gd_k = []
        for (i, j) in zip(gK, dK):
            gd_k.append(i * j)
        sum_gd_K = sum(gd_k)

        #Segunda parte de la desigualdad
        secondPart = 0.5 * ak * sum_gd_K
        if firstPart <= secondPart:
            break

        else:
            ak /= 2
    print(ak)
    return ak

# Calcula bk
def calc_bk(gk, pre_gk, dk):
    
    # Se calcula la norma 2 del vector actual
    gk_norma = linalg.norm(array(gk, dtype='float'), 2)

    # Se calcula la norma 2 del vector anterior
    pre_gk_norma = linalg.norm(array(pre_gk, dtype='float'), 2)

    bk = gk_norma ** 2 / pre_gk_norma ** 2
    return bk
 

gradiente("(x-2)**4 + (x-2*y)**2", ['x', 'y'], [0, 3], 10**(-10), 13) 