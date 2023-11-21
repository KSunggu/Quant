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


xw.Books()