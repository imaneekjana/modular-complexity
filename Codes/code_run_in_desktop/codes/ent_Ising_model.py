#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 22:14:55 2024

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

L=12 # system size
J=1.0 # spin zz interaction

def create_ent(h):
    
    
    h_field=[[-h,i] for i in range(L)]
    J_zz=[[-J,i,(i+1)%L] for i in range(L)] # PBC

    static_spin =[["zz",J_zz],["x",h_field]] # static part of H
    dynamic_spin=[] # time-dependent part of H

    basis = spin_basis_1d(L=L) 

    H_spin=hamiltonian(static_spin,dynamic_spin,basis=basis,dtype=np.float64)

    E,V = H_spin.eigh()

    state = np.outer( V.T[0], conj(V.T[0].T))

    ent_spec = basis.ent_entropy(state,sub_sys_A= range(basis.L//2),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                                sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)["p_A"]
    
    return ent_spec


#h_list = np.arange(0.1,2,0.1) 

h_list = [0.25,0.5,1.0,1.25,1.5,0.75]
tf = 200
dt = 0.1
t_list = np.arange(0,tf,dt)


ent_list = []

complexity_R_diff_h =[]
complexity_L_diff_h =[]

capacity_ent = []

for hval in h_list:
    lamb = create_ent(hval)
    modular_E = [- np.log(i) for i in lamb ]
    lancz = lanczos_modular(modular_E)
    ent = lancz['an_list'][0]
    CE = (lancz['bn_list'][0])**2
    ent_list.append(ent)
    capacity_ent.append(CE)
    
    
    comp_list =  modular_complexity(modular_E,tf,dt)
    complR = comp_list['compR']
    complL = comp_list['compL']
    np.savetxt("/home/user/Documents/Official/Courses/krylov/krylov_moduler/TFIM_modular/L12/R_modular_complexity_J_{}_h_{}_L_{}.txt".format(J,hval,L),complR)
    np.savetxt("/home/user/Documents/Official/Courses/krylov/krylov_moduler/TFIM_modular/L12/L_modular_complexity_J_{}_h_{}_L_{}.txt".format(J,hval,L),complL)
    complexity_R_diff_h.append(complR)
    complexity_L_diff_h.append(complL)
    
    
np.savetxt("/home/user/Documents/Official/Courses/krylov/krylov_moduler/TFIM_modular/L12/t_list.txt",t_list)

'''
np.savetxt("/home/user/Documents/Official/Courses/krylov/krylov_moduler/TFIM_modular/L12/h_list.txt",h_list)  

np.savetxt("/home/user/Documents/Official/Courses/krylov/krylov_moduler/TFIM_modular/L12/entanglement_entropy_gr_state_J_{}_L_{}.txt".format(J,L),ent_list)  

np.savetxt("/home/user/Documents/Official/Courses/krylov/krylov_moduler/TFIM_modular/L12/capacity_of_entanglement_gr_state_J_{}_L_{}.txt".format(J,L),capacity_ent)  
'''      
    
    
  
plt.figure()
plt.plot(h_list,ent_list,linestyle='-.',label = 'ent',color='b')
plt.plot(h_list,capacity_ent,linestyle='-.',label = 'CE',color='k')
plt.xlabel('h')
plt.legend()
plt.show()   

'''    
plt.figure()
plt.plot(t_list,complexity_R_diff_h[0],label = 'compR, h = 0.25',color='r')
#plt.plot(t_list,complexity_L_diff_h[0],label = 'compL, h = 0.25',color='k')
plt.plot(t_list,complexity_R_diff_h[1],label = 'compR, h = 0.5',color='b')
plt.plot(t_list,complexity_R_diff_h[2],label = 'compR, h = 1.0',color='k')
plt.plot(t_list,complexity_R_diff_h[3],label = 'compR, h = 1.25',color='g')
#plt.plot(t_list,complexity_L_diff_h[1],label = 'compL, h = 1.0',color='k')
plt.plot(t_list,complexity_R_diff_h[4],label = 'compR, h = 1.75',color='y')
plt.plot(t_list,complexity_R_diff_h[5],label = 'compR, h = 0.75',color='c')
#plt.plot(t_list,complexity_L_diff_h[2],label = 'compL, h = 1.75',color='k')
plt.legend()
plt.show()    '''
    













