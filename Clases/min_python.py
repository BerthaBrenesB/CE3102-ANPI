# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 16:30:56 2021

@author: Dell
"""

import sympy as sp
from scipy import optimize

f='x^2+x-5'
x=sp.Symbol('x')
fs=sp.sympify(f)
fn=sp.lambdify(x,fs)
num_min = optimize.fmin(fn, 1)
"""Acá se obtiene el valor de x donde se alcanza el mínimo"""
print(num_min[0]) 