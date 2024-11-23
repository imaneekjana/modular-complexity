#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:34:17 2024

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

np.conjugate = conj

L=8

basis = spin_basis_1d(L=L) 

## The function below, returns the ground state given the value of J, and h

def create_gs(J,Jp,h,g):
    
    #L=10 # system size
    J_list = [J if i != (L//2 -1) else Jp for i in range (L-1)]
    J_zz=[[J_list[i],i,(i+1)] for i in range(L-1)] # OBC
    h_field=[[h,i] for i in range(L)]
    g_field=[[g,i] for i in range(L)]
    

    static_spin =[["zz",J_zz],["z",h_field],["x",g_field]] # static part of H
    dynamic_spin=[] # time-dependent part of H

    H_spin=hamiltonian(static_spin,dynamic_spin,basis=basis,dtype=np.float64)

    E,V = H_spin.eigh()
    
    return V.T[0]


J1 = -1
Jp1 = 0
h1 = 0.1
g1 = 1.5

J2 = -1
Jp2 = -1
h2 = 0.1
g2 = 1.5

J3 = -1
Jp3 = 1
h3 = 0.1
g3 = 1.5


gs1 = create_gs(J1,Jp1,h1,g1) # before quencing the Hamiltonians
gs2 = create_gs(J2,Jp2,h2,g2) # after local quench
gs3 = create_gs(J3,Jp3,h3,g3)

#print(gs1)

rho_f12 = np.outer(gs1,gs2.T.conj())
rho_f11 = np.outer(gs1,gs1.T.conj())
rho_f22 = np.outer(gs2,gs2.T.conj())
rho_f33 = np.outer(gs3,gs3.T.conj())
rho_f13 = np.outer(gs1,gs3.T.conj())


rho_A12 = (basis.ent_entropy(rho_f12,sub_sys_A= np.arange(L//4,3*L//4,1),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                           sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A12 = rho_A12/np.trace(rho_A12)
lam_spec12 = np.linalg.eigvals(rho_A12)
e_spec12o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec12]

e_spec12 = []

for en in e_spec12o:
    if en !=0:
        e_spec12.append(en)
        
        
rho_A13 = (basis.ent_entropy(rho_f13,sub_sys_A= np.arange(L//4,3*L//4,1),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                           sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A13 = rho_A13/np.trace(rho_A13)
lam_spec13 = np.linalg.eigvals(rho_A13)
e_spec13o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec13]

e_spec13 = []

for en in e_spec13o:
    if en !=0:
        e_spec13.append(en)        

#ent_12 = sum([(-lam*np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec12])


rho_A11 = (basis.ent_entropy(rho_f11,sub_sys_A= np.arange(L//4,3*L//4,1),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                            sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A11 = rho_A11/np.trace(rho_A11)
lam_spec11 = np.linalg.eigvals(rho_A11)
e_spec11o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec11]

e_spec11 = []

for en in e_spec11o:
    if en !=0:
        e_spec11.append(en)

#ent_11 = sum([(-lam*np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec11])


rho_A22 = (basis.ent_entropy(rho_f22,sub_sys_A= np.arange(L//4,3*L//4,1),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                            sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A22 = rho_A22/np.trace(rho_A22)
lam_spec22 = np.linalg.eigvals(rho_A22)
e_spec22o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec22]

e_spec22 = []

for en in e_spec22o:
    if en !=0:
        e_spec22.append(en)
        
        
        
rho_A33 = (basis.ent_entropy(rho_f33,sub_sys_A= np.arange(L//4,3*L//4,1),return_rdm="A",enforce_pure=False,return_rdm_EVs=True,density = False,
                            sparse=False,alpha=1.0,sparse_diag=True,subsys_ordering=True)['rdm_A'])
rho_A33 = rho_A33/np.trace(rho_A33)
lam_spec33 = np.linalg.eigvals(rho_A33)
e_spec33o = [(-np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec33]

e_spec33 = []

for en in e_spec33o:
    if en !=0:
        e_spec33.append(en)        
#ent_22 = sum([(-lam*np.log(lam) if abs(lam) > 10**(-15) else 0) for lam in lam_spec22])

#print(np.real(ent_12)-0.5*(ent_11+ent_22))

os.chdir("/home/user/Documents/Official/Courses/krylov/krylov_moduler/TFIM_with_longitudinal_field/L8")
os.makedirs("Pseudo_mod_comp_g_{}_h_{}_J_{}".format(g1,h1,J1))

tf = 500
dt = 0.1

t_list = np.arange(0,tf,dt)

np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/t_list.txt".format(g1,h1,J1),t_list)

comp12 = modular_complexity(e_spec12, tf, dt)
comp12_R = comp12['compR']
comp12_L = comp12['compL']

comp13 = modular_complexity(e_spec13, tf, dt)
comp13_R = comp13['compR']
comp13_L = comp13['compL']

comp11 = modular_complexity(e_spec11, tf, dt)
comp11_R = comp11['compR']
comp11_L = comp11['compL']

comp22 = modular_complexity(e_spec22, tf, dt)
comp22_R = comp22['compR']
comp22_L = comp22['compL']


comp33 = modular_complexity(e_spec33, tf, dt)
comp33_R = comp33['compR']
comp33_L = comp33['compL']





np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp11_R.txt".format(g1,h1,J1),comp11_R)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp11_L.txt".format(g1,h1,J1),comp11_L)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp22_R.txt".format(g1,h1,J1),comp22_R)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp22_L.txt".format(g1,h1,J1),comp22_L)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp33_R.txt".format(g1,h1,J1),comp33_R)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp33_L.txt".format(g1,h1,J1),comp33_L)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp12_R.txt".format(g1,h1,J1),comp12_R)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp12_L.txt".format(g1,h1,J1),comp12_L)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp13_R.txt".format(g1,h1,J1),comp13_R)
np.savetxt("Pseudo_mod_comp_g_{}_h_{}_J_{}/comp13_L.txt".format(g1,h1,J1),comp13_L)

plt.figure()

plt.plot(t_list,comp12_R,label='compR12')
plt.plot(t_list,comp12_L,label='compL12')

plt.plot(t_list,comp11_R,label='compR11')
#plt.plot(t_list,comp11_L,label='compL11')

plt.plot(t_list,comp22_R,label='compR22')
#plt.plot(t_list,comp22_L,label='compL22')
plt.plot(t_list,comp33_R,label='compR33')
plt.plot(t_list,comp13_R,label='compR13')
plt.legend()

plt.title(f'g_{g1}_h_{h1}')

plt.show()


