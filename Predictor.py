#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 01:22:13 2020

@author: bs15ansj
"""

from sympy import sympify
import numpy as np


class Model:
    def __init__(self, a=10, polymer_model='gaussian_chain',
                 restriction_factor=2*np.pi, linker_receptor_dG=0):
        self.a = a
        self.polymer_model = 'gaussian_chain'
        self.restriction_factor = restriction_factor
        self.kds = []

    
    def calculate_kds(self, kd, rmin=10, rmax=40, step=1, micro_kd=0):
    
        for r in range(rmin, rmax+1):

            if kd == 'macro':
                answers = self.macro_polynomial(r)
                for answer in answers:
                    if sympify(answer).is_real and answer > 0.0:
                        self.kds.append(answer)
                        print(answer)
            elif kd == 'micro':
                self.kds.append(self.micro_polynomials(r, micro_kd))
                
                    
    def save_kds(self, fname):
        with open(fname, 'w+') as f:
            for i, kd in enumerate(self.kds):
                f.write(str(kd)+'\n')
                print(str(i)+' done; ' + str(kd))
    
    
            
            
