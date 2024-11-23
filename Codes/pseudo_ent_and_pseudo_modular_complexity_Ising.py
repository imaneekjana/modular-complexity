#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 18:44:04 2024

@author: maitri
"""

from __future__ import print_function, division
from quspin.operators import hamiltonian, exp_op  # Hamiltonians and operators
from collections import Counter
from scipy.sparse import csr_matrix  # for sparse matrix tranformation
import matplotlib.pyplot as plt  # plotting library
import numpy as np  # generic math functions
from quspin.tools.block_tools import block_diag_hamiltonian  # block diagonalisation
from quspin.basis import spinless_fermion_basis_1d  # Hilbert space fermion basis
import sys
import os
qspin_path = os.path.join(os.getcwd(), "../../")
sys.path.insert(0, qspin_path)
import cmath
from scipy.stats import gaussian_kde

from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis import spin_basis_1d, spinless_fermion_basis_1d # Hilbert space spin basis
import numpy as np # generic math functions
import matplotlib.pyplot as plt # plotting library

from krylov_modular import *

##### define model parameters #####
np.conjugate = conj

L=10

basis = spin_basis_1d(L=L) 

## The function below, returns the ground state given the value of J, and h

def create_gs(J,h):
    
    #L=10 # system size
    
    
    h_field=[[-h,i] for i in range(L)]
    J_zz=[[-J,i,(i+1)%L] for i in range(L)] # PBC

    static_spin =[["zz",J_zz],["x",h_field]] # static part of H
    dynamic_spin=[] # time-dependent part of H

    H_spin=hamiltonian(static_spin,dynamic_spin,basis=basis,dtype=np.float64)

    E,V = H_spin.eigh()
    
    return V.T[0]


## finding reduced density matrix

J1= 1
h1 = 1.51

J2 = 1
h2 = 0.51

gs1 = create_gs(J1, h1)
gs2 = create_gs(J2,h2)

#print(gs1)

rho_f12 = np.outer(gs1,gs2.T.conj())
rho_f11 = np.outer(gs1,gs1.T.conj())
rho_f22 = np.outer(gs2,gs2.T.conj())


rho_A12 = (basis.ent_entropy(rho_f12,sub_sys_A= range(basis.L//2),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                           sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A12 = rho_A12/np.trace(rho_A12)
lam_spec12 = np.linalg.eigvals(rho_A12)
e_spec12o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec12]

e_spec12 = []

for en in e_spec12o:
    if en !=0:
        e_spec12.append(en)

#ent_12 = sum([(-lam*np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec12])


rho_A11 = (basis.ent_entropy(rho_f11,sub_sys_A= range(basis.L//2),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                            sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A11 = rho_A11/np.trace(rho_A11)
lam_spec11 = np.linalg.eigvals(rho_A11)
e_spec11o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec11]

e_spec11 = []

for en in e_spec11o:
    if en !=0:
        e_spec11.append(en)

#ent_11 = sum([(-lam*np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec11])


rho_A22 = (basis.ent_entropy(rho_f22,sub_sys_A= range(basis.L//2),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                            sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A22 = rho_A22/np.trace(rho_A22)
lam_spec22 = np.linalg.eigvals(rho_A22)
e_spec22o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec22]

e_spec22 = []

for en in e_spec22o:
    if en !=0:
        e_spec22.append(en)
#ent_22 = sum([(-lam*np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec22])

#print(np.real(ent_12)-0.5*(ent_11+ent_22))


tf = 20
dt = 0.1

t_list = np.arange(0,tf,dt)

comp12 = modular_complexity(e_spec12, tf, dt)
comp12_R = comp12['compR']
comp12_L = comp12['compL']

comp11 = modular_complexity(e_spec11, tf, dt)
comp11_R = comp11['compR']
comp11_L = comp11['compL']

comp22 = modular_complexity(e_spec22, tf, dt)
comp22_R = comp22['compR']
comp22_L = comp22['compL']



plt.figure()

plt.plot(t_list,comp12_R,label='compR12')
plt.plot(t_list,comp12_L,label='compL12')

plt.plot(t_list,comp11_R,label='compR11')
plt.plot(t_list,comp11_L,label='compL11')

plt.plot(t_list,comp22_R,label='compR22')
plt.plot(t_list,comp22_L,label='compL22')

plt.legend()

plt.title(f'$J_1={J1},h_1={h1}\,\,\,J_2={J2},h_2={h2}$')

plt.show()


















