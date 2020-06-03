#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 00:57:31 2020

@author: bs15ansj
"""


import Systems
import numpy as np
import matplotlib.pyplot as plt

'''
This script uses predicts the Kd of hetero-bivalent inhibitors of CTB as a 
function of linker length for flexible and rigid linkers and plots the results
'''

def calculator():
    '''
    This function calculates Kd as a function of linker length for flexible 
    and bivalent ligands binding to GM1 and LeY of CTB
    '''
    model = Systems.CTB(a=5, polymer_model='gaussian_chain')
    model.calculate_kds('macro', rmin=10, rmax=40, step=0.01)
    model.save_kds('example/flex.dat')
    
    model = Systems.CTB(a=5, polymer_model='harmonic_spring')
    model.calculate_kds('macro', rmin=10, rmax=40, step=0.01)
    model.save_kds('example/rigid.dat')
    

def plotter():
    
    '''
    This function plots the results from calculator()
    '''
    
    r = np.linspace(10.0,40.0,31)
    
    flex = np.loadtxt('example/flex.dat')
    rigid = np.loadtxt('example/rigid.dat')
    
    plt.plot(r, flex, color='red', linestyle='-', label='flex')
    plt.plot(r, rigid, color='blue', linestyle='-', label='rigid')
    plt.plot(r, r-r+43.3e-9, color='black', linestyle=':', linewidth=4, label='GM1')

    plt.legend()
    plt.ylim(10e-12, 10e-6)
    plt.xlim(10,40)
    plt.ylabel(r'$K_{D}$' + ' (M)', rotation='horizontal')
    plt.xlabel(r'$r_{ete}$' + r' ($\AA$)')
    plt.grid(axis='y')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.yscale('log')
    plt.savefig('example/plot.png', dpi=300)

calculator()
plotter()