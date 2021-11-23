import math
import sympy as sp 
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

def Lk(xv, k):
    x = sp.Symbol('x')
    [n,m] = xv.shape
    Lk = 1
    for i in range(0,n):
        if(i != k):
            Lk = Lk*( (x -xv[i]/(xv[k]-xv[i])))
    return Lk

        

def lagrange(xk,yk):
    xk = np.matrix(xk)
    yk = np.matrix(yk)
    x = sp.Symbol('x')
    [n,m] = xk.shape
    print(n)
    poly = 0;
    for k in range(0,n):
        poly += yk[k] * Lk(xk,k)
    return poly

xk = '-2; 0; 1';
xk_5 = '0; 2; 4;6;8';
yk = '0; 1; -1';
yk_5 = '-0.2; -0.6; -1.4;-3;-6.2';
p = lagrange(xk_5, yk_5);
print(p)
        