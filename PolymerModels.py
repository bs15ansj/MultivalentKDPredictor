#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 01:13:46 2020

@author: bs15ansj
"""

import numpy as np

def Ceff_Gaussian(a, r, d):
    t1a = 8.0 * np.sqrt(6.0) * (a**2)
    t1b = (np.pi**(3.5)) * (r**5)
    t1 = t1a/t1b
    t2a = -1.5
    t2b = (d/r)**2
    t2 = np.exp(t2a*t2b)
    ans = (t1*t2)*(10000/6.022)
    return ans
        
def Ceff_Harmonic(d, r, delta_r):
    t1 = 1.0/(2.0*(np.pi)**(3.0/2.0))
    t2 = np.exp((-1.0/2.0)*(((d-r)**2.0)/(delta_r**2.0)))
    t3 = delta_r*r**2.0
    ans = (t1*(t2/t3))*(10000.0/6.022)
    return ans