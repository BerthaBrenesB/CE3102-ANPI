# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 16:19:13 2021

@author: Usuario
"""

import numpy as np

A=np.array([[-3,1,0],[1,-2,1],[0,1,-3]])
"npA=np.linalg.eig(A)"
"print(npA)"
x0=np.array([[1],[1],[1]])
iterMax=100
tol=10**-10

xk=x0

for k in range(0,iterMax):
    yk=np.dot(A,xk)
    ck=np.linalg.norm(yk,np.inf)
    xk_n=(1/ck)*yk
    error=np.linalg.norm(xk_n-xk)
    xk=xk_n 
    if error<tol:
        break
    
print('Valor Propio Mayor =', ck)
print('Vector propio respectivo =', xk)
print('Iteraciones =', k)
print('Error =', error)