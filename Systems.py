#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 01:00:27 2020

@author: bs15ansj
"""

import numpy as np
from sympy import Sum, symbols, binomial, sympify
from sympy.solvers import solve

from Predictor import Model
import PolymerModels

D, i, ig, il, ia, ib = symbols('D i ig il ia ib')

class CTB(Model):

    def __init__(self, a, polymer_model, restriction_factor=2*np.pi, linker_receptor_dG=0):
        super().__init__(a, polymer_model, restriction_factor=2*np.pi, linker_receptor_dG=0)
        self.gm1_kd = 43.3e-9 # https://pubs.acs.org/doi/10.1021/ja0378207
        self.ley_kd = 1.4e-3 # https://pubs.acs.org/doi/10.1021/acsinfecdis.7b00085
        self.gm1_ka = 1/self.gm1_kd
        self.ley_ka = 1/self.ley_kd
        self.polymer_model = polymer_model 
        self.dleya = 24.1 # from pymol
        self.dleyb = 32.8 # from pymol
        self.a = a
        self.restriction_factor = restriction_factor
        self.linker_receptor = np.exp(linker_receptor_dG)

    
    def macro_polynomial(self, r):
        if self.polymer_model == 'gaussian_chain':
            leya_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Gaussian(self.a, r, self.dleya)*self.restriction_factor*self.linker_receptor
            leyb_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Gaussian(self.a, r, self.dleyb)*self.restriction_factor*self.linker_receptor
        elif self.polymer_model =='harmonic_spring':
            leya_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Harmonic(self.dleya, r, r/10)*self.restriction_factor*self.linker_receptor
            leyb_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Harmonic(self.dleyb, r, r/10)*self.restriction_factor*self.linker_receptor
        roots = solve(self.gm1(self.gm1_ka, r) + self.leya(leya_ka)
                      + self.leyb(leyb_ka) + self.gm1_and_leya(leya_ka)
                      + self.gm1_and_leyb(leyb_ka)  
                      + self.leya_and_leyb(leya_ka, leyb_ka)
                      + self.gm1_and_leya_and_leyb(leya_ka, leyb_ka)
                      - 5, D)
        
        return roots
    
    def micro_polynomials(self, r, micro_kd):
        if self.polymer_model == 'gaussian_chain':
            leya_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Gaussian(self.a, r, self.dleya)*self.restriction_factor*self.linker_receptor
            leyb_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Gaussian(self.a, r, self.dleyb)*self.restriction_factor*self.linker_receptor
        elif self.polymer_model =='harmonic_spring':
            leya_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Harmonic(self.dleya, r, r/10)*self.restriction_factor*self.linker_receptor
            leyb_ka = self.gm1_ka*self.ley_ka*PolymerModels.Ceff_Harmonic(self.dleyb, r, r/10)*self.restriction_factor*self.linker_receptor
        
        if micro_kd == 0:
            return 1/leya_ka
        elif micro_kd == 1:
            return 1/leyb_ka

    def gm1(self, Kgm1, r):
        if self.polymer_model == 'gaussian_chain':
            alpha = 6*((2/(3*np.pi))**(3/2))*(self.a/r)
        elif self.polymer_model == 'harmonic_spring':
            alpha = 0.5
        else:
            print('no alpha!')
            
        Kgm1 = Kgm1*alpha*self.linker_receptor
        gm1 = (2*i-5)*binomial(5, i)*(Kgm1**i*D**i)
        gm1 = Sum(gm1, (i, 1, 5)).doit()
        return gm1

    def leya(self, Kleya):
        leya = (2*i-5)*binomial(5, i)*(Kleya**i*D**i)
        leya = Sum(leya, (i, 1, 5)).doit()
        return(leya)
    
    def leyb(self, Kleyb):
        leyb = (2*i-5)*binomial(5, i)*(Kleyb**i*D**i)
        leyb = Sum(leyb, (i, 1, 5)).doit()
        return(leyb)
    
    def gm1_and_leya(self, Kleya):
        gm1_and_leya = ((2*ig) + (2*ia) - 5)*binomial(5, ig)*binomial(5-ig, ia)*(self.gm1_ka**ig)*(Kleya**ia)*D**(ig+ia)
        gm1_and_leya = Sum(Sum(gm1_and_leya, (ig, 1, 5)).doit(), (ia, 1, 5)).doit()
        return(gm1_and_leya)
    
    def gm1_and_leyb(self, Kleyb):
        gm1_and_leyb = ((2*ig) + (2*ib) - 5)*binomial(5, ig)*binomial(5-ig, ib)*(self.gm1_ka**ig)*(Kleyb**ib)*D**(ig+ib)
        gm1_and_leyb = Sum(Sum(gm1_and_leyb, (ig, 1, 5)).doit(), (ib, 1, 5)).doit()
        return(gm1_and_leyb)
    
    def leya_and_leyb(self, Kleya, Kleyb):
        leya_and_leyb1 = ((2*ia)+2-5)*binomial(5,1)*binomial(3,ia)*Kleya**ia*Kleyb*D**(ia+1)
        leya1_and_leyb = ((2*ib)+2-5)*binomial(5,1)*binomial(3,ib)*Kleya*Kleyb**ib*D**(ib+1)
        leya2_and_leyb2 = (4+4-5)*binomial(5,1)*Kleya**2*Kleyb**2*D**(4)
        leya_and_leyb = Sum(leya_and_leyb1, (ia, 1, 3)).doit() + Sum(leya1_and_leyb, (ib, 1, 3)).doit() + leya2_and_leyb2
        return leya_and_leyb
    
    def gm1_and_leya_and_leyb(self, Kleya, Kleyb):
        gm1_and_leya_and_leyb1 = ((2*ia)+2+(2*ig)-5)*binomial(5,1)*binomial(3,ia)*binomial(4-ia,ig)*self.gm1_ka**ig*Kleya**ia*Kleyb*D**(ig+ia+1)
        gm1_and_leya1_and_leyb = ((2*ib)+2+(2*ig)-5)*binomial(5,1)*binomial(3,ib)*binomial(4-ib,ig)*self.gm1_ka**ig*Kleya*Kleyb**ib*D**(ig+ib+1)
        gm1_leya2_and_leyb2 = (4+4+2-5)*binomial(5,1)*self.gm1_ka*Kleya**2*Kleyb**2*D**(5)
        gm1_leya_and_leyb = (Sum(Sum(gm1_and_leya_and_leyb1, (ia, 1, 3)).doit(), (ig, 1, 5)).doit()
                          + Sum(Sum(gm1_and_leya1_and_leyb, (ib, 1, 3)).doit(), (ig, 1, 5)).doit()
                          + gm1_leya2_and_leyb2)
        return gm1_leya_and_leyb
    
    


    
    