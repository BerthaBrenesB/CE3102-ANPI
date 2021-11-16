function archivo_predictor_corrector
  clc;
  pkg load symbolic;
  syms x y; #Variable a trabajar 
  warning('off', 'all');
  f =  'y -x^2 + 1'; #Ecuacion diferencial para ser resuelta por el metodo
  funcion = matlabFunction(sym(f)); #Se convierte a simbolico la funcion
  intervalo = [0,2]; #Intervalo donde se define la funcion 
  y0 = 0.5; #Valor inicial en y
  puntos = 11; #Cantidad de puntos 
  predictor_corrector(funcion,y0,intervalo,puntos); #Se llama la funcion predictor_corrector
end


#         FUNCION PREDICTOR CORRECTOR
#
#     Parametros de Entrada:
#     funcion = Funcion para ser resuelta por el metodo
#     y0 = Valor inicial de y para la funcion
#     intervalo = Intervalo donde se define la funcion
#     puntos = Puntos donde trabaja el metodo para la funcion
#
#     Parametros de Salida:
#     polinomio = El polinomio de interpolacion con los x y y de la funcion obtenidos por el metodo
#     Grafica de puntos 
function polinomio = predictor_corrector(funcion,y0, intervalo, puntos)
  
  a = intervalo(1); #Se define el valor de a 
  b = intervalo(2); #Se define el valor de b
  
  h = ( b - a ) / (puntos - 1); #Se define el valor de h para el metodo
  xn = a:h:b; #Se incializa el valor de x para las iteraciones 
  yn = [y0]; #Se inicializa el valor de y para la iteraciones
  
  for n=1:puntos-1
    
    #Seccion predictor
    y_n_p = yn(n)+(h*funcion(xn(n), yn(n))); #Se procede al calculo de yn+1 con Euler
    numerador = funcion(xn(n), yn(n)) + funcion(xn(n +1), y_n_p); #S resuelve el yn+1
    #Seccion corrector
    y_n_c = yn(n) + ( h * (numerador / 2)); #Calculo del yn+1 predictor corrector
    
    yn = [yn y_n_c]; #Se almacenan los valores de y en las itreaciones del metodo
    
  end
  
    polinomio = metodo_lagrange(xn, yn) #Se calcula el polinomio de interpolacion
    
    plot(xn,yn,'b')  #Se grafica el comportamiento de la funcion 
    title("Grafica de puntos de comportamiento") #Se coloca el titulo a la grafica 
    xlabel("Valores en x (XV)") #Se nombre el eje x 
    ylabel("Valores en y (YV)") #Se nombra el eje y 
end


#                 FUNCION LK
#
#     Parametros de Entrada:
#     xv = Se refiere a la lista de valores en x 
#     k = Numero de la iteracion donde trabaja el metodo
#
#     Parametros de Salida:
#     LK = Valores de los factores para el calculo de Lagrange
#     Grafica de puntos 
function Lk=fun_Lk(xv,k)
  syms x  #Se define el simbolico para la funcion
  n=length(xv)-1;  #Se toma el largo de la lista menos uno
  Lk=1;   #Se incializa el valore de los factores en 1 
  for j=0:n   #Se incializa el ciclo
    if j~=k   #Si son diferentes 
      Lk=Lk*(x-xv(j+1))/(xv(k+1)-xv(j+1));  #Se aplica el calculo del factor 
    end    
  end
  Lk=expand(Lk); #Se obtiene el resultado
end

#                 FUNCION METODO LAGRANGE 
#
#     Parametros de Entrada:
#     xv = Se refiere a la lista de valores en x 
#     YV = Se refiere a la lista de valores en y
#
#     Parametros de Salida:
#     p = Polimonio de interpolacion asociado a los valores previos en x y y
function p=metodo_lagrange(xv,yv)
  syms x    #Se define el simbolico para la funcion
  n=length(xv)-1; #Se toma el largo de la lista menos uno 
  p=0;  #Se inicializa el polinomio en menos uno 
  for k=0:n #Se incializa el ciclo 
    p=p+yv(k+1)*fun_Lk(xv,k);  #Se procede con el calculo del polinomio de interpolacion
  end
  p=expand(p);  #Se obtiene el resultado del polinomio 
end



   
   
