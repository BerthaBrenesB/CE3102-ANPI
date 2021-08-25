%Ejemplo de evaluar un numero a en la 
%función exponencial, utilizando el polinomio 
%de Taylor
clc; clear;

tol=1e-15;
a=2;
Sk=0;
iterMax=100000;
e=[];
for n=0:iterMax
  Sk_n=Sk+a^n/factorial(n);
  error=abs(Sk_n-Sk)/abs(Sk_n);
  e=[e error];
  if error<tol
    break
  end
  Sk=Sk_n;
end

plot(0:length(e)-1,e,'g','LineWidth',2)
title('Grado del Polinomio vrs Error Relativo')
xlabel('Grado del Polinomio (k)')
ylabel('Error Relativo')
