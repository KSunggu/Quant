# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 21:47:59 2023

@author: Sunggu
"""

import QuantLib as ql
import numpy as np
import xlwings as xw
import matplotlib.pyplot as plt

# HW 1factor model dr(t) = a(theta(t)/a - r(t))dt + sigma dz
# theta(t) = F_t(0,t) + aF(0,t) + sigma^2

sigma = 0.01995;


file = xw.books("back data.xlsx")
Inputs = file.sheets('Sheet1')


t = Inputs.range('A2').expand('down').value
r = Inputs.range('B2').expand('down').value

a = 10;
vol = 0.01;
r0 = 0.03;

plt.plot(t, r, 'o-')
plt.title("Yield Curve")
plt.grid(visible = True)




def make_spline(x, a):
    n = len(x) -1;
    h = []
    alpha = np.zeros(n);
    for i in range(n):
        h.append(x[i+1] - x[i]);
        
    for i in range(1, n):
        alpha[i] = 3/ h[i] *(a[i+1] - a[i]) -3 /h[i-1] *(a[i] - a[i-1]);
        
    l = np.zeros(n+1);
    l[0] = 1;
    z =np.zeros(n+1)
    for i in range(1,n):
        l[i] = 2*(x[i+1] - x[i-1]) -h[i-1] - 
    