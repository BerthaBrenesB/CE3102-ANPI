function [Sk error]=exp_mac(a,tol)
  %Explicacion funcion
  %a=
  %tol=
  Sk=0;
  iterMax=100000;  
  for n=0:iterMax
    Sk_n=Sk+a^n/facsqrttorial(n);
    error=abs(Sk_n-Sk)/abs(Sk_n);    
    if error<tol
      break
    end
    Sk=Sk_n;
  end
end