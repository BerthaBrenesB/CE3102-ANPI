function calculoTest()
  f1 = '(x-2)^4 + (x-2*y)^2'
  x0 = [0, 3]
  tol = 10^-10
  iterMax=2
  [xk_1] = GNCL(f1,x0,iterMax,tol)
end;

function [xk_1, iter]=GNCL(f, x0, iterMax, tol)
%    Funcion que implementa el metodo de gradiente conjugado
%    f: string con la funcion que se debe evaluar, debe tener solamente dos variables
%    x0: vector inicial que necesita el metodo
%    tol: tolerancia al fallo que debe tener el resultado
%    retorna una lista de dos elementos, el calculo de xk+1 y numero de iteraciones
  pkg load symbolic
% Se obtitne las variables y la funcion simbolica
  syms x y;
  f = sym(f);
% Se obtiene una version numerica de la funcion para ser evaluada
  f1 = matlabFunction(f);
% Se obtiene el gradiente simbolico de la funcion
  g = gradient(f);
% se obtiene la version numerica de la gradiente 
  g_n = matlabFunction(g);
% El valor de gk, el cual es el gradiente de f evaluado con el vector inicial 
  gk = g_n(x0(1),x0(2))'
  dk = -gk
  xk=x0
% Variables for the plot
  I = []
  E = []
  for k=0:iterMax
 % calculo de ak
    ak = calculo_ak1(f1,xk,dk,gk);
 % Calculo de xk+1 = xk + ak*dk
    xk_1 = xk + ak.*dk
    I(end+1) = k;
    E(end+1) = norm(g_n(xk_1(1),xk_1(2)));
% Calculo del error el cual es la norma del gradiente de f(x+1)
    if(norm(g_n(xk_1(1),xk_1(2)))<tol);
      xk = xk_1;
     break; % si el valor del error es menor a la tolerancia
    endif
 % una vez calculado el error se pueden calcular nuevamente las gk+1 y bk
    gk_1 = g_n(xk_1(1), xk_1(2))';
    bk = calculo_bk(gk_1,gk);
 % El nuevo calculo de dk+1
    dk_1 = -gk_1 + bk.*dk;
  endfor
  plot(I,E)
  
endfunction

function [ak] = calculo_ak1(f,xk,dk,gk)
%    Funcion que calcula el valor de ak
%    f: funcion simbolica que se debe evaluar, debe tener solamente dos variables
%    xk: vector inicial que necesita el metodo
%    dk: el valor del vector dk
%    gk; el valor del vector gk
%    retorna el valor de ak
% Primero se calcula con un ak =1 y l = 0.5
  ak = 1;
  while (true)
    % calculo del valor a evaluar la funcion en la parte izq f(xk + ak*dk)
    xk_adk = xk + ak*dk;
    % Calculo de la parte izq de la desigualdad
    izq = f(xk_adk(1),xk_adk(2)) -f(xk(1),xk(2));
    % Calculo de la parte derecha de la desigualdad
    der = 0.5*ak* sum(dk.*gk);
    % Evaluacion de la igual si se cumple se retorna el ak 
    if(izq <= der)
      break;
    else
    % Si no, puede seguir evaluando con ak menor al anterior
      ak = ak/2;
    endif
  endwhile
endfunction

function [bk] = calculo_bk(gk, pre_gk)
%    Funcion que calcula el valor de bk
%    gk; el valor del vector gk+1
%    pre_gk: el valor de gk
%    retorna el valor de bk
  bk = norm(gk)^2/norm(pre_gk)^2;
endfunction

