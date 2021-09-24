clc; clear; 
%A=[1 2 3 4; 0 1 2 1; 2 4 6 8;0 -3 -6 -3]

m=30;
A=rand(m,m/2)*rand(m/2,m);

for k=1:m  
  Ak=A(1:k,1:k);  
  deter=det(Ak)
  tol=eps;
  if abs(deter)<tol
    k    
    break
  endif
end
