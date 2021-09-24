# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 15:02:48 2021

@author: Dell
"""

import sympy as sp
import numpy as np

def variables_simbolicas(variables):
    n=len(variables)
    tam=np.arange(0,n,2)
    v_var=[]
    for i in tam:
        v_var.append(sp.Symbol(variables[i]))
    return v_var
        
def gradiente(f,v_var):
    n=len(v_var)
    g=[]
    for i in np.arange(0,n,1):
        g.append(sp.diff(f,v_var[i]))
    return g


"""Ejemplo de Variables"""
variables='x y'
v_var=variables_simbolicas(variables)
print(v_var)

"""Ejemplo del Gradiente Simbolico"""
f='(x-2)^4+(x-2*y)^2'
f_s=sp.sympify(f)
g=gradiente(f_s,v_var)
print(g)

"""Ejemplo de Sustitución Simbolica y Numérica"""
x0=0.65
y0=0.25
f_n=sp.lambdify(v_var,f)
print(f_n(x0,y0))
print(f_s.subs([(v_var[0],x0),(v_var[1],y0)]))

g_n=sp.lambdify(v_var,g)
print(g_n(x0,y0))