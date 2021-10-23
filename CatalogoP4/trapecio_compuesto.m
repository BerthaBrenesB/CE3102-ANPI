function archivo_trapecio_compuesto
  pkg load symbolic
  f='log(x)' 
  intervalo=[2,5]
  num = 500
  [aprox,error]=trapecio_compuesto(f,num,intervalo)
  
end
function [aprox,error]=trapecio_compuesto(f,n,intervalo)
  f = sym(f)
  f1 = matlabFunction(f)
  a = intervalo(1);
  b = intervalo(2);
  h = (b-a)/(n-1)
  x0 = a;
  xv = linspace(a,b,n);
  I=0;
  for i=1:n-1
    ai = xv(i);
    bi = xv(i+1);
    fai = f1(ai,'x');
    fbi = f1(bi,'x');
    I += ((bi-ai)*(fai+fbi))/2;
  endfor
  aprox = I;
  error = cota_error_trapecio(f,intervalo,h)
endfunction
function [error]=cota_error_trapecio(f, intervalo,h)
  a = intervalo(1);
  b = intervalo(2);
  f2d = abs(diff(diff(f,'x')))
  fn2d = matlabFunction(f2d)
  d_fa = fn2d(a,'x');
  d_fb = fn2d(b,'x');
  d2_fx = 0
  if(d_fa> d_fb)
    d2_fx = d_fa;
  else
    d2_fx = d_fb;
  end
  error = (((b-a)*h**2/12))*d2_fx;
endfunction
