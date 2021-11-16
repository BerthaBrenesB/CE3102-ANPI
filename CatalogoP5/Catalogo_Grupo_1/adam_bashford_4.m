function archivo_adam_bashford_4
  clc; 
  pkg load symbolic;
  syms x y;
  #pkg load symbolic;
  warning('off', 'all');
  funcion =  '1 + (x - y)^2'; #Ecuacion diferencial para ser resuelta por el metodo
  intervalo = [2,4]; #Intervalo donde se define la funcion 
  y0 = 1; #Valor inicial en y
  x0 = 2; #Valor inicial en x
  y1 = 1.191; #Valor y_1 en y 
  puntos = 11; #Cantidad de puntos 
  adam_bashford_4(funcion,y0, x0, y1, intervalo,puntos); #Se llama la funcion predictor_corrector
end


#         FUNCION ADAM BASHFORD 4
#
#     Parametros de Entrada:
#     funcion = Funcion para ser resuelta por el metodo
#     y0 = Valor inicial de y para la funcion
#     x0 = Valor inicial de x para la funcion
#     y1 = Valor siguiente para y1
#     intervalo = Intervalo donde se define la funcion
#     puntos = Puntos donde trabaja el metodo para la funcion
#
#     Parametros de Salida:
#     polinomio = El polinomio de interpolacion con los x y y de la funcion obtenidos por el metodo
#     Grafica de puntos 
function polinomio = adam_bashford_4(funcion,y0, x0, y1, intervalo,puntos)
  a = intervalo(1); #Se define el valor de a 
  b = intervalo(2); #Se define el valor de b
  
  h = ( b - a ) / (puntos - 1); #Se define el valor de h para el metodo
  
  listaX = [x0, (x0 + h), (x0 + 2*h), (x0 + 3*h) ]; #Valores inciales para x
  listaY = [y0, y1, 1.365, 1.528]; #Valores inciales para y
  
  for contador=4:puntos-1
  
    #Se obtiene el valor de y_temp
    listaY(contador+1) =listaY(contador)+(h/24)*(55*funcion(listaX(contador),listaY(contador)) - 59*funcion(listaX(contador-1),listaY(contador-1)) + 37*funcion(listaX(contador-2),listaY(contador-2)) - 9*funcion(listaX(contador-3),listaY(contador-3)));
    listaX(contador+1) = listaX(contador) + h #Se incrementa el valor de x_temp en h 

  end
  
   polinomio = metodo_lagrange(listaX, listaY) #Se calcula el polinomio de interpolacion
    
   plot(xn,yn,'b')  #Se grafica el comportamiento de la funcion 
   title("Grafica de puntos de comportamiento") #Se coloca el titulo a la grafica 
   xlabel("Valores en x (ListaX)") #Se nombre el eje x 
   ylabel("Valores en y (ListaY)") #Se nombra el eje y 

  
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