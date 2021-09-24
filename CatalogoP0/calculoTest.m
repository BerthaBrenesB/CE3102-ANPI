function calculoTest()
  f1 = '(x-2)^4 + (x-2*y)^2';
  %xk = [2.05551301351897, 1.02768000760174];
  %dk = [-0.00132822742613825, 0.000891155474249500];
  %gk= [0.000990293262482567, -0.000611993261967925];
  x0 = [0, 3]
  tol = 10^-10
  iterMax=2
  [xk_1] = GNCL(f1,x0,iterMax,tol)
end;

function [xk_1, iter]=GNCL(f, x0, iterMax, tol)
  pkg load symbolic
  syms x y;
  f = sym(f);
  f1 = matlabFunction(f);
  var = ['x', 'y']
  g = gradient(f);
  g_n = matlabFunction(g);
  gk = g_n(x0(1),x0(2))'
  dk = -gk
  xk=x0
  iter = 0
  I = []
  E = []
  for k=0:iterMax
    ak = calculo_ak1(f1,xk,dk,gk);
    xk_1 = xk + ak.*dk;
    I(end+1) = k;
    E(end+1) = norm(g_n(xk_1(1),xk_1(2)));
    if(norm(g_n(xk_1(1),xk_1(2)))<tol);
      xk = xk_1;
      break;
    endif
    gk_1 = g_n(xk_1(1), xk_1(2))';
    bk = calculo_bk(gk_1,gk);
    dk_1 = -gk_1 + bk.*dk;
    iter = iter+1;
  endfor
  plot(I,E)
  
endfunction

function [ak] = calculo_ak1(f,xk,dk,gk)
  ak = 1;
  while (true)
    xk_adk = xk + ak*dk;
    izq = f(xk_adk(1),xk_adk(2)) -f(xk(1),xk(2));
    der = 0.5*ak* sum(dk.*gk);
    if(izq <= der)
      break;
    else
      ak = ak/2
    endif
  endwhile
endfunction

function [bk] = calculo_bk(gk, pre_gk)
  bk = norm(gk)^2/norm(pre_gk);
endfunction

