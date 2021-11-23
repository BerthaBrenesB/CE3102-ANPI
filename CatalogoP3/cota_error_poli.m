% Cota de Error Polinomio de Interpolacion
clc; clear;
pkg load symbolic;
syms x;

function [cotaError] = cota_poly_inter(f, s)
%entrada:
%funcion y puntos
%salida:
%cota error

a = s(1); %valor minimo intervalo
b = s(end); %valor maximo de intervalo
val = 0.54; %punto a evaluar

n = length(s) - 1;
multValor = 1;
% formula para evaluar el punto val - resta
  for k=0:n
    restValor = val - s(k+1);
    multValor = multValor * restValor;
  endfor
  absolutoMult = abs(multValor);
  
  fs=sym(f);
  n_1 = n + 1;
  derivada_n1 = diff(fs,n_1)
  faux = -1*abs(derivada_n1)
  fauxNum = matlabFunction(faux); %convierte a funcion numerica 
  
  xMax = fminbnd(fauxNum, a, b) %valor en x donde la funcion derivada se hace maximo
  cotaError = (xMax * absolutoMult) / factorial(n_1) %calculo de la cota
  
end
% prueba del metodo para la funcion
%p = (4*x)/3 -(x^3)/3; %polinomio
f = sin((pi*x)/2)
s = [-1, 0, 1, 2];
[cotaError] = cota_poly_inter(f, s);