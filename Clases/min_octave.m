%Ejemplo minimo de una funcion en una variable

clc; clear;
pkg load symbolic

f='x^2-3*x+20';
fs=sym(f);
fn=matlabFunction(fs);
x0=1;
[xmin,ymin]=fminsearch(fn,x0)