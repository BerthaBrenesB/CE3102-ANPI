function ejemplo_biseccion()
 % Calcular el cero de exp(x) -2*x-10 =0
 % PAs 1: Conocer el intervalo donde se encuentra el cero
 % sugerencia: Graficar la funcion 
 %f(x) = 2^-x^
 %$xv=-5:5;
 %yv= exp(xv) -2*xv -10;
 %plot(xv,yv);
 %grid on
 
 % Concluimos que la funcion tiene un 0 en el intervalo [2,4] a partir de la grafica
 %Cadeba de caracteres
 %f = 'exp(x) -2*x-10 ';
 f= '2 + cos(e**x -2) -e**x';
 %Extreo a y b
 a = 1;
 b = 9;
 %tol un numero positivo que representa la toletacion
 tol = 10^-5;
 % Cantidad de iteraciones maximas
 iterMax = 10;
 %Retorno con los variables de salida
 [xk k error] = biseccion(f,a,b,tol, iterMax)
 
end;

function [xk k error] = biseccion(f,a,b,tol,iterMax)
  %Esta función aproxima la solución de la ecuación f(x)=0, 
  % utilizando el método de la bisección
  %
  %Sintaxis:  [xk k error]=biseccion(f,a,b,tol,iterMax)
  % 
  %Parámetros Iniciales: 
  %            f = una  cadena de caracteres (string) que representa a la función f
  %            a,b = son los extremos del intervalo [a,b]
  %            tol = un número positivo que representa a la tolerancia para el criterio |f(xk)|<tol
  %            iterMax = cantidad de iteraciones máximas
  %            
  %Parámetros de Salida:                           
  %            xk = aproximación del cero de la función f
  %            k = número de iteraciones realizados
  %            error =  |f(xk)|
  
  %%%% Se debe instalar el paquete Symbolic
  %%%% Paso 1: Descargar el archivo symbolic-win-py-bundle-2.9.0.tar.gz
  %%%%         de la pagina https://github.com/cbm755/octsympy/releases
  %%%% Paso 2: Escribir en la Venta de Comandos de Octave la instruccion
  %%%%         pkg install symbolic-win-py-bundle-2.9.0.tar.gz  
  
  %%% Cargar el paquete symbolic
  pkg load symbolic
  f1 = matlabFunction(sym(f));
  if f1(x,a)*f1(x,b)<0
    %se cumple el teorema de bolzano
    for k=1:iterMax
      xk=(a+b)/2; %Tenemos 2 intervales: [a xk] y [xk b]
      if f1(x,a)*f1(x,xk)<0 %Se cumple la condicion en el intervalo 1
        b=xk;
      else % Se cumple la condicion en el intervalo 2
        a=xk; 
      end
      error= abs(f1(xk));
      if error<tol
        break
      end
    end
  else
    xk='NA';
    k='NA';
    error='NA';
    display('El intervalo selecionado no cumple las condiciones del teorema de Bolzano');
  end
end
