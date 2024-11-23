#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 16:39:13 2024

@author: maitri
"""

from krylov_modular import *

import numpy as np
 
import matplotlib.pyplot as plt

import scipy as sp

import cmath

from scipy.special import *

np.conjugate = conj


def SU2_lambda(z1,z2,l):
    lambda_list = []
    for n in range(int(2*l+1)):
        lamb = (1 + z1*conj(z2))**(-2*l) *( gamma(2*l+1)/(gamma(2*l-n+1)*factorial(n)))*(z1*conj(z2))**n
        lambda_list.append(lamb)
    return lambda_list    
        

### pseduo-entropy excess
'''
delS_list = []
delC_list = []
phi_list = np.arange(-7,7,np.pi/50)

for phi in phi_list:
    
    z1 = 1.0
    z2 = 1.0*cmath.exp(1j*phi)
    
    E_list12 = [-np.log(lam) for lam in SU2_lambda(z1, z2, 1/2)]
    
    E_list11 = [-np.log(lam) for lam in SU2_lambda(z1, z1, 1/2)]
    
    E_list22 = [-np.log(lam) for lam in SU2_lambda(z2, z2, 1/2)]
    
    S12 = lanczos_modular(E_list12)['an_list'][0]
    S11 = lanczos_modular(E_list11)['an_list'][0]
    S22 = lanczos_modular(E_list22)['an_list'][0]
    
    C12 = (lanczos_modular(E_list12)['bn_list'][0])**2
    C11 = (lanczos_modular(E_list11)['bn_list'][0])**2
    C22 = (lanczos_modular(E_list22)['bn_list'][0])**2
    
    
    
    delS = np.real(S12)- 0.5*(S11+S22)
    
    delC = np.real(C12)- 0.5*(C11+C22)
    
    
    delS_list.append(delS)
    delC_list.append(delC)
    
plt.figure()

plt.plot(phi_list,delS_list,label='$\Delta S$')

plt.plot(phi_list,delC_list,label='$\Delta C$')

plt.legend()
plt.show()

'''
### modular complexity

z1 = 0.4

z2 = 1.1+0.4j

tf = 50.0
dt=0.01

t_list = np.arange(0,tf,dt)

E_list12 = [-np.log(lam) for lam in SU2_lambda(z1, z2, 5/2)]

comp = modular_complexity(E_list12, tf, dt)
compR = comp['compR']
compL = comp['compL']

plt.figure()

plt.plot(t_list,compR,label='compR')
plt.plot(t_list,compL,label='compL')

plt.legend()

plt.show()
    

