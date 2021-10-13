function Lk=fun_Lk(xv,k)
  syms x
  %k=0,1,....,n
  n=length(xv)-1;
  Lk=1;
  for j=0:n
    if j~=k
      Lk=Lk*(x-xv(j+1))/(xv(k+1)-xv(j+1))
  end
  Lk=expand(Lk);
end
