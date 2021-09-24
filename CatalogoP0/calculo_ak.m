function calculoTest()
  f1 = '(x-2)^4 + (x-2*y)^2';
  %xk = [2.05551301351897, 1.02768000760174];
  %dk = [-0.00132822742613825, 0.000891155474249500];
  %gk= [0.000990293262482567, -0.000611993261967925];
  x0 = [0, 3]
  tol = 10^-10
  iterMax=13
  [xk] = GNCL(f1,x0,iterMax,tol)
end;

function [xk]=GNCL(f, x0, iterMax, tol)
  pkg load symbolic
  f = matlabFunction(sym(f));
  g0 = gradient(f(x0))
  d0 = -g0
  xk=x0
  for k=0:iterMax
    ak = calculo_ak1(f, xk,dk,gk)
    x(k+1) = xk + ak.*dk
    if(norm(gradient(f(xk))) < tol)
      xk = x(k+1)
      break
    endif
    g(k+1) = gradient(f(xk+1))
    bk = calculo_bk(g(k+1), gk)
    d(k+1) = -gk+1 + bk.*dk
  endfor
endfunction
function [ak] = calculo_ak1(f,xk,dk,gk)
  ak = 1
  while (true)
    xk_adk = xk + ak*dk
    izq = f(xk_adk(1),xk_adk(2)) -f(xk(1),xk(2))
    der = 0.5*ak* sum(dk.*gk)
    if(izq <= der)
      break;
    else
      ak = ak/2
    endif
  endwhile
endfunction

function [bk] = calculo_bk(gk, gk_1)
  bk = norm(gk)^2/norm(gk_1)
endfunction
