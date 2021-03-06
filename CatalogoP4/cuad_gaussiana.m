function archivo_cuad_gaussiana
  clc;
  pkg load symbolic
  warning('off', 'all');

  funcion='exp(x)*cos(x)'; %Ejemplo de como definir los datos
  intervalo=[-2,2]; %Intervalo donde se define la integral 
  orden=4; %Orden de la funcion para aplicar el metodo
  [error,aprox]=cuad_gaussiana(funcion,orden,intervalo)#Resultado de la aproximacion
end


%                   FUNCION CUADRATURA GAUSSIANA 
%     Funcion encargada de calcular el valor de una integral por el metodo de cuadratura gaussiana
%
%     Parametros de entrada
%     funcion:  funcion a evaluar para la integral
%     invervalo:  Intervalo donde se define la integral
%     orden: orden de derivada para la funcion a aproximar
%
%     Parametros de salida
%     aproximacion: valor de la aproximacion con el metodo trapecio, 
%     error: valor de la cota de error
function [error,aprox]=cuad_gaussiana(funcion,orden,intervalo)
  funcion=sym(funcion); %Se transforma la funcion de string a simbolica
  x=sym('x'); %Se determina a x como la variable
  a=intervalo(1); %Valor minimo del intervalo 
  b=intervalo(2); %Valor maximo del intervalo 
  g_x=((b-a)/2)*subs(funcion,((b-a)*x+(b+a))/2); %Formula para la ecuacion para la aproximacion en el intervalo dado
  aprox=cuad_gaussiana_aux(g_x,orden);#Se calcula la aproximacion con respecto a g
  e1=cuad_gaussiana_aux(funcion,orden);#Se calcula la aproximacion integral de -1 a 1 de f
  error=abs(aprox-e1);#Calculo de valor de error
end 

%                   FUNCION CUADRATURA GAUSSIANA AUXILIAR
%     Funcion encargada de calcular el valor de una integral por el metodo de cuadratura gaussiana
%     El objetivo de esta funcion es poder acceder a los ceros y pesos de las derivadas de la funcion
%     segun sea su orden, esto para evitar el calculo de la derivada
%
%     Parametros de entrada
%     funcion:  funcion a evaluar para la integral
%     invervalo:  Intervalo donde se define la integral
%     orden: orden de derivada para la funcion a aproximar
%
%     Parametros de salida
%     aproximacion: valor de la aproximacion con el metodo trapecio, 
%     error: valor de la cota de error
function [aproximacion]=cuad_gaussiana_aux(funcion,orden)
  funcion=matlabFunction(funcion);#Convierte la funcion a tipo matlab
  %x corresponde a los ceros segun el orden 
  %w corresponde a los pesos segun el orden 
  switch orden#Valores de la derivada ya definidos segun el orden
    case 2
      x=[ -0.577350269189626  0.577350269189626];
      w=[1 1];
    case 3
      x=[ -0.774596669241483 0  0.774596669241483];
      w=[ 0.555555555555556 0.888888888888889 0.555555555555556];
    case 4  
      x=[-0.86113631159405 -0.339981043584856 0.339981043584856 0.86113631159405];
      w=[0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454];
    case 5
      x=[-0.906179845938664 -0.538469310105683 0 0.538469310105683 0.906179845938664];
      w=[0.236926885056189 0.478628670499366 0.568888888888889 0.478628670499366 0.236926885056189];
    case 6
      x=[-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
      w=[0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    case 7
      x=[-0.949107912342759 -0.741531185599394 -0.405845151377397 0 0.405845151377397 0.741531185599394 0.949107912342759];
      w=[0.129484966168870 0.279705391489277 0.381830050505119 0.417959183673469 0.381830050505119 0.279705391489277 0.129484966168870];
    case 8
      x=[-0.960289856497536 -0.796666477413627 -0.525532409916329 -0.183434642495650 0.183434642495650 0.525532409916329 0.796666477413627 0.960289856497536];
      w=[0.101228536290376 0.222381034453374 0.313706645877887 0.362683783378362 0.362683783378362 0.313706645877887 0.222381034453374 0.101228536290376];
    case 9
      x=[-0.968160239507626 -0.836031107326636 -0.613371432700590 -0.324253423403809 0 ...
         0.324253423403809 0.613371432700590 0.836031107326636 0.968160239507626];
      w=[0.081274388361574 0.180648160694857 0.260610696402935 0.312347077040003 0.330239355001260 ...
         0.312347077040003 0.260610696402935 0.180648160694857 0.081274388361574];
    case 10
      x=[-0.973906528517172 -0.865063366688985 -0.679409568299024 -0.433395394129247 -0.148874338981631 ...
         0.148874338981631 0.433395394129247 0.679409568299024 0.865063366688985 0.973906528517172];
      w=[0.066671344308688 0.149451349150581 0.219086362515982 0.269266719309996 0.295524224714753 ... 
          0.295524224714753 0.269266719309996 0.219086362515982 0.149451349150581 0.066671344308688];
    otherwise
      error("El mayor orden corresponde a 10")
  endswitch
  aproximacion=0; %Se inicializa la variable para el calculo del error
  #Se itera sobre el orden dado para el calculo del metodo
  for i=1:orden
    aproximacion=aproximacion+w(i)*funcion(x(i));#Calcula la aproximacion por el metodo de cuadratura gaussiana 
  endfor
end