xv = [-2 0 1]
yv = [0 1 -1]
p = metodo_lagrange(xv,yv)

function p=metodo_lagrange(xv,yv)
  syms x
  n=length(xv)-1;
  p=0;
  for k=0:n
    p=p+yv(k+1)*fun_Lk(xv,k);
  end
  p=expand(p);
end
