#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 20:55:46 2024

@author: maitri
"""

import numpy as np

import matplotlib.pyplot as plt

import scipy as sp

conj = np.conjugate
#---------------------------------------Creating the action of H and H^+ where H is the Hamiltonian---------------------------------------------


def actH(state,E_data):
    dim = len(E_data)
    state_f = np.zeros(dim, dtype = np.complex128)
    for i in range(dim):
        st_ith = E_data[i]*state[i]
        state_f[i] = st_ith
        
    return state_f

def actH_dag(state,E_data):
    dim = len(E_data)
    state_f = np.zeros(dim, dtype = np.complex128)
    for i in range(dim):
        st_ith = conj(E_data[i])*state[i]
        state_f[i] = st_ith
        
    return state_f   
        
def actU(state,E_data,dt):
    dim = len(E_data)
    state_f = np.zeros(dim, dtype = np.complex128)
    
    for i in range(dim):
        st_ith = np.exp(-1j*E_data[i]*dt)*state[i]
        state_f[i] = st_ith
        
    return state_f

def actU_dag(state,E_data,dt):
    dim = len(E_data)
    state_f = np.zeros(dim, dtype = np.complex128)
    
    for i in range(dim):
        st_ith = np.exp(-1j*conj(E_data[i])*dt)*state[i]
        state_f[i] = st_ith
        
    return state_f 
#-------------------------------------Creating the an,bn,cn,Pn,Qn lists-----------------------------------------------------------------------------------

def lanczos_modular(E_data):
    dim = len(E_data)
    P0 = np.zeros(dim, dtype = np.complex128)
    for i in range(dim):
        P0[i] = np.exp(-E_data[i]/2)
        
    Q0 = np.zeros(dim, dtype = np.complex128)    
    for i in range(dim):
        Q0[i] = np.exp(-conj(E_data[i])/2)    
    a0 = conj(Q0).T @ actH(P0,E_data)
    
    an_list = [a0]
    bn_list = []
    cn_list = []
    Pn_list = [P0]
    Qn_list = [Q0]
    
    for n in np.arange(1,dim,1):
        if n==1:
            A1 = actH(Pn_list[-1], E_data) - an_list[-1]*Pn_list[-1]
            B1 = actH_dag(Qn_list[-1], E_data) - conj(an_list[-1])*Qn_list[-1]
            w1 = conj(B1).T@A1
            b1 = np.sqrt(w1)
            c1 = b1
            P1 = A1/c1
            Q1 = B1/conj(b1)
            a1 = conj(Q1).T @ actH(P1,E_data)
            if np.abs(b1) > 10**(-6):
                an_list.append(a1)
                bn_list.append(b1)
                cn_list.append(c1)
                Pn_list.append(P1)
                Qn_list.append(Q1)
            else:
                break
            
                
            
            
            
        else: 
            An = actH(Pn_list[-1], E_data) - an_list[-1]*Pn_list[-1] - bn_list[-1]*Pn_list[-2]
            Bn = actH_dag(Qn_list[-1], E_data) - conj(an_list[-1])*Qn_list[-1] - conj(cn_list[-1])*Qn_list[-2]
            
            proj_A = np.zeros(dim,dtype = np.complex128)
            proj_B = np.zeros(dim,dtype = np.complex128)
            for m in range(n):
                proj_A += (conj(Qn_list[m]).T @ An)*Pn_list[m]
                proj_B += (conj(Pn_list[m]).T @ Bn)*Qn_list[m]
            
            An = An - proj_A 
            Bn = Bn - proj_B 
            
            proj_A = np.zeros(dim,dtype = np.complex128)
            proj_B = np.zeros(dim,dtype = np.complex128)
            for m in range(n):
                proj_A += (conj(Qn_list[m]).T @ An)*Pn_list[m]
                proj_B += (conj(Pn_list[m]).T @ Bn)*Qn_list[m]
            
            An = An - proj_A 
            Bn = Bn - proj_B 
            
            wn = conj(Bn).T@An
            bn = np.sqrt(wn)
            cn = bn
            Pn = An/cn
            Qn = Bn/conj(bn)
            an = conj(Qn).T @ actH(Pn,E_data)
            if np.abs(bn) > 10**(-6):
                an_list.append(an)
                bn_list.append(bn)
                cn_list.append(cn)
                Pn_list.append(Pn)
                Qn_list.append(Qn)
            else:
                break
            
            
            
    Dict =  {} 
    Dict['an_list'] = an_list
    Dict['bn_list'] = bn_list
    Dict['cn_list'] = cn_list
    Dict['Pn_list'] = Pn_list
    Dict['Qn_list'] = Qn_list
    
    return Dict
        
#------------------------------------- Calculate the complexity ------------------------------------------------

def modular_complexity(E_data,tf,dt):
    t_list =  np.arange(0,tf,dt)
    
    Dict = lanczos_modular(E_data)
    
    Pn_list = Dict['Pn_list']
    Qn_list = Dict['Qn_list']
    
    Psi_R = Pn_list[0]
    Psi_L = Qn_list[0]
    
    PsiR_working = Psi_R
    PsiL_working = Psi_L
    complexity_R_list = []
    complexity_L_list = []
    for t in t_list:
        PsiR_t = PsiR_working/(np.linalg.norm(PsiR_working))
        PsiL_t = PsiL_working/(np.linalg.norm(PsiL_working))
        
        complexity_R = 0
        complexity_L = 0
        norm_R = 0
        norm_L = 0
        
        for n in range(len(Pn_list)):
            phi_R = conj(Qn_list[n]).T @ PsiR_t
            phi_L = conj(Pn_list[n]).T @ PsiL_t
            
            norm_L += np.abs(phi_L)**2
            norm_R += np.abs(phi_R)**2
            
            complexity_L += n*np.abs(phi_L)**2
            complexity_R += n*np.abs(phi_R)**2
            
            
        complexity_L = complexity_L/norm_L  
        complexity_R = complexity_R/norm_R 
        
        complexity_R_list.append(complexity_R)
        complexity_L_list.append(complexity_L)
        
        PsiR_working = actU(PsiR_t,E_data,dt)
        PsiL_working = actU_dag(PsiL_t,E_data,dt)
    
    Dict_comp = {}
    Dict_comp['compR'] = np.array(complexity_R_list)
    Dict_comp['compL'] = np.array(complexity_L_list)
    
    return Dict_comp

    


#--------------------------------------------------------------------------------------------
tf = 100
t_list = np.arange(0,tf,0.01)
dt = 0.01
p = 1/3 + 1j*0.05  

q = 1/3 + 1j*0.05

E_data = [-np.log(p),-np.log(q), - np.log(1-p-q)]


modon = modular_complexity(E_data, tf, dt)

plt.figure()
plt.plot(t_list,modon['compR'],label='R')   
plt.plot(t_list,modon['compL'],label='L')    
plt.legend()    
plt.show()        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        