#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 21:00:05 2024

@author: maitri
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex']=True
'''
CE_10 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L10/capacity_of_entanglement_gr_state_J_1.0_L_10.txt",dtype = np.complex128)
SE_10 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L10/entanglement_entropy_gr_state_J_1.0_L_10.txt",dtype = np.complex128)
CE_12 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/capacity_of_entanglement_gr_state_J_1.0_L_12.txt",dtype = np.complex128)
SE_12 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/entanglement_entropy_gr_state_J_1.0_L_12.txt",dtype = np.complex128)

h_val = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L10/h_list.txt",dtype = np.complex128)

plt.figure()
plt.plot(h_val,CE_10,label = "$ CE, L = 10 $",linestyle='-.')
plt.plot(h_val,CE_12,label = "$ CE, L = 12 $",linestyle='-.')
plt.plot(h_val,SE_10,label = "$ SE, L = 10 $",linestyle='-.')
plt.plot(h_val,SE_12,label = "$ SE, L = 12 $",linestyle='-.')
plt.xlabel("$ h $", fontsize=22)
#plt.ylabel("$ CE $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=4.0)
plt.xticks(np.arange(0,2.1,0.5),fontsize=22) 
plt.yticks(np.arange(0,1,0.2),fontsize=22)   
plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/CE_SE_gr_state_J_1.0_L_10.pdf")
'''
#---------------------------------------------------------------------------------------------------------------------------------------------

L_mc_1 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/L_modular_complexity_J_1.0_h_0.25_L_12.txt",dtype = np.complex128)
L_mc_2 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/L_modular_complexity_J_1.0_h_0.5_L_12.txt",dtype = np.complex128)
L_mc_3 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/L_modular_complexity_J_1.0_h_0.75_L_12.txt",dtype = np.complex128)
L_mc_4 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/L_modular_complexity_J_1.0_h_1.0_L_12.txt",dtype = np.complex128)
L_mc_5 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/L_modular_complexity_J_1.0_h_1.25_L_12.txt",dtype = np.complex128)
L_mc_6 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/L_modular_complexity_J_1.0_h_1.5_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/t_list.txt")

plt.figure()
plt.plot(t_list,L_mc_1,label = "$ h = 0.25 $")
plt.plot(t_list,L_mc_2,label = "$ h = 0.50 $")
plt.plot(t_list,L_mc_3,label = "$ h = 0.75 $")
plt.plot(t_list,L_mc_4,label = "$ h = 1.00 $")
plt.plot(t_list,L_mc_5,label = "$ h = 1.25 $")
plt.plot(t_list,L_mc_6,label = "$ h = 1.50 $")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_L  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(loc='upper left', bbox_to_anchor=(1, 1),fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
#plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=3.0)
plt.xticks(np.arange(0,201,50),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/left_modular_complexity_L12.pdf")

R_mc_1 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/R_modular_complexity_J_1.0_h_0.25_L_12.txt",dtype = np.complex128)
R_mc_2 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/R_modular_complexity_J_1.0_h_0.5_L_12.txt",dtype = np.complex128)
R_mc_3 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/R_modular_complexity_J_1.0_h_0.75_L_12.txt",dtype = np.complex128)
R_mc_4 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/R_modular_complexity_J_1.0_h_1.0_L_12.txt",dtype = np.complex128)
R_mc_5 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/R_modular_complexity_J_1.0_h_1.25_L_12.txt",dtype = np.complex128)
R_mc_6 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/R_modular_complexity_J_1.0_h_1.5_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/t_list.txt")

plt.figure()
plt.plot(t_list,R_mc_1,label = "$ h = 0.25 $")
plt.plot(t_list,R_mc_2,label = "$ h = 0.50 $")
plt.plot(t_list,R_mc_3,label = "$ h = 0.75 $")
plt.plot(t_list,R_mc_4,label = "$ h = 1.00 $")
plt.plot(t_list,R_mc_5,label = "$ h = 1.25 $")
plt.plot(t_list,R_mc_6,label = "$ h = 1.50 $")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_R  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(loc='upper left', bbox_to_anchor=(1, 1),fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
#plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=3.0)
plt.xticks(np.arange(0,201,50),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/right_modular_complexity_L12.pdf")

#-----------------------------------------------------------------------------------------------------------------------------------------------
'''
L_pmc_11_a = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity11_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt",dtype = np.complex128)
L_pmc_11_b = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity11_J1_1_h1_0.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)
L_pmc_11_c = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity11_J1_1_h1_0.25_J2_1_h2_0.5_L_12.txt",dtype = np.complex128)
L_pmc_11_d = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity11_J1_1_h1_1.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/t_list_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt")

plt.figure()
plt.plot(t_list, L_pmc_11_a, label="$ h_1 = 0.5, h_2 = 1.5$")
plt.plot(t_list, L_pmc_11_b, label="$ h_1 = 0.5, h_2 = 1.75$")
plt.plot(t_list, L_pmc_11_c, label="$ h_1 = 0.25, h_2 = 0.5$")
plt.plot(t_list, L_pmc_11_d, label="$ h_1 = 1.5, h_2 = 1.75$")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_L^{11}  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=4.0)
plt.xlim(0,800)
plt.xticks(np.arange(0,800,150),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()

plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/left_pseudo_modular_complexity11_L12.pdf")
#---------------------------------------------------------------------------------------------------------------------------------------------------


L_pmc_22_a = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity22_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt",dtype = np.complex128)
L_pmc_22_b = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity22_J1_1_h1_0.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)
L_pmc_22_c = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity22_J1_1_h1_0.25_J2_1_h2_0.5_L_12.txt",dtype = np.complex128)
L_pmc_22_d = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity22_J1_1_h1_1.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/t_list_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt")

plt.figure()
plt.plot(t_list, L_pmc_22_a, label="$ h_1 = 0.5, h_2 = 1.5$")
plt.plot(t_list, L_pmc_22_b, label="$ h_1 = 0.5, h_2 = 1.75$")
plt.plot(t_list, L_pmc_22_c, label="$ h_1 = 0.25, h_2 = 0.5$")
plt.plot(t_list, L_pmc_22_d, label="$ h_1 = 1.5, h_2 = 1.75$")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_L^{22}  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=4.0)
plt.xlim(0,800)
plt.xticks(np.arange(0,800,150),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/left_pseudo_modular_complexity22_L12.pdf")
#------------------------------------------------------------------------------------------------------------------------------------------------

L_pmc_12_a = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity12_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt",dtype = np.complex128)
L_pmc_12_b = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity12_J1_1_h1_0.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)
L_pmc_12_c = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity12_J1_1_h1_0.25_J2_1_h2_0.5_L_12.txt",dtype = np.complex128)
L_pmc_12_d = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/L_pseudo_modular_complexity12_J1_1_h1_1.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/t_list_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt")

plt.figure()
plt.plot(t_list, L_pmc_12_a, label="$ h_1 = 0.5, h_2 = 1.5$")
plt.plot(t_list, L_pmc_12_b, label="$ h_1 = 0.5, h_2 = 1.75$")
plt.plot(t_list, L_pmc_12_c, label="$ h_1 = 0.25, h_2 = 0.5$")
plt.plot(t_list, L_pmc_12_d, label="$ h_1 = 1.5, h_2 = 1.75$")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_L^{12}  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=4.0)
plt.xlim(0,800)
plt.xticks(np.arange(0,800,150),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()

plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/left_pseudo_modular_complexity12_L12.pdf")'''

#---------------------------------------------------------------------------------------------------------------------------------------------
'''
R_pmc_11_a = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity11_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt",dtype = np.complex128)
R_pmc_11_b = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity11_J1_1_h1_0.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)
R_pmc_11_c = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity11_J1_1_h1_0.25_J2_1_h2_0.5_L_12.txt",dtype = np.complex128)
R_pmc_11_d = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity11_J1_1_h1_1.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/t_list_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt")

plt.figure()
plt.plot(t_list, R_pmc_11_a, label="$ h_1 = 0.5, h_2 = 1.5$")
plt.plot(t_list, R_pmc_11_b, label="$ h_1 = 0.5, h_2 = 1.75$")
plt.plot(t_list, R_pmc_11_c, label="$ h_1 = 0.25, h_2 = 0.5$")
plt.plot(t_list, R_pmc_11_d, label="$ h_1 = 1.5, h_2 = 1.75$")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_R^{11}  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=4.0)
plt.xlim(0,800)
plt.xticks(np.arange(0,800,150),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()

plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/right_pseudo_modular_complexity11_L12.pdf")
#---------------------------------------------------------------------------------------------------------------------------------------------------


R_pmc_22_a = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity22_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt",dtype = np.complex128)
R_pmc_22_b = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity22_J1_1_h1_0.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)
R_pmc_22_c = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity22_J1_1_h1_0.25_J2_1_h2_0.5_L_12.txt",dtype = np.complex128)
R_pmc_22_d = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity22_J1_1_h1_1.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/t_list_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt")

plt.figure()
plt.plot(t_list, R_pmc_22_a, label="$ h_1 = 0.5, h_2 = 1.5$")
plt.plot(t_list, R_pmc_22_b, label="$ h_1 = 0.5, h_2 = 1.75$")
plt.plot(t_list, R_pmc_22_c, label="$ h_1 = 0.25, h_2 = 0.5$")
plt.plot(t_list, R_pmc_22_d, label="$ h_1 = 1.5, h_2 = 1.75$")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_R^{22}  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=4.0)
plt.xlim(0,800)
plt.xticks(np.arange(0,800,150),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/right_pseudo_modular_complexity22_L12.pdf")
#------------------------------------------------------------------------------------------------------------------------------------------------

R_pmc_12_a = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity12_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt",dtype = np.complex128)
R_pmc_12_b = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity12_J1_1_h1_0.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)
R_pmc_12_c = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity12_J1_1_h1_0.25_J2_1_h2_0.5_L_12.txt",dtype = np.complex128)
R_pmc_12_d = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/R_pseudo_modular_complexity12_J1_1_h1_1.5_J2_1_h2_1.75_L_12.txt",dtype = np.complex128)

t_list = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/t_list_J1_1_h1_0.5_J2_1_h2_1.5_L_12.txt")

plt.figure()
plt.plot(t_list, R_pmc_12_a, label="$ h_1 = 0.5, h_2 = 1.5$")
plt.plot(t_list, R_pmc_12_b, label="$ h_1 = 0.5, h_2 = 1.75$")
plt.plot(t_list, R_pmc_12_c, label="$ h_1 = 0.25, h_2 = 0.5$")
plt.plot(t_list, R_pmc_12_d, label="$ h_1 = 1.5, h_2 = 1.75$")

plt.xlabel("$ s $", fontsize=22)
plt.ylabel("$ \mathcal C_R^{12}  $", fontsize=22)
# Add borders to the plot
ax = plt.gca()
# Set the border color and width
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Set the border color
    spine.set_linewidth(2.5)      # Set the border width
    
plt.legend(fontsize = 18, frameon = False, fancybox = False, shadow = False, borderpad = 0)
plt.gca().set_aspect('auto', adjustable='box')
plt.tight_layout(pad=4.0)
plt.xlim(0,800)
plt.xticks(np.arange(0,800,150),fontsize=22) 
plt.yticks(np.arange(0,3,0.5),fontsize=22)   
plt.show()

plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/right_pseudo_modular_complexity12_L12.pdf")'''

#--------------------------------------------------------------------------------------------------------------------------------------
'''
Delta_CE_pseudo_L12 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/Delta_CE_J1_1_J2_1_L_12.txt",dtype = np.complex128)

Delta_SE_pseudo_L12 = np.loadtxt("/home/maitri/Documents/krylov_moduler/TFIM_modular/L12/pseudo/Delta_SE_J1_1_J2_1_L_12.txt",dtype = np.complex128)


plt.figure()
im = plt.imshow(np.real(Delta_CE_pseudo_L12), cmap = "RdBu",interpolation="nearest",origin='lower')
#plt.colorbar(label="$\Delta C_E$",fontsize= 22)
plt.xticks(np.arange(0,101,20),fontsize=22) 
plt.yticks(np.arange(0,101,20),fontsize=22)  
# Add a colorbar and set font sizes
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=20)  # Change font size of colorbar ticks
cbar.set_label("$\Delta C_E$", fontsize=22)  # Change font size of the colorbar label

plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/Delta_CE_pseudo_L12.pdf")

plt.figure()
im = plt.imshow(abs(np.real(Delta_CE_pseudo_L12)), cmap = "magma",interpolation="nearest",origin='lower')
#plt.colorbar(label="$\Delta C_E$",fontsize= 22)
plt.xticks(np.arange(0,101,20),fontsize=22) 
plt.yticks(np.arange(0,101,20),fontsize=22)  
# Add a colorbar and set font sizes
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=20)  # Change font size of colorbar ticks
cbar.set_label("$|\Delta C_E|$", fontsize=22)  # Change font size of the colorbar label

plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/abs_Delta_CE_pseudo_L12.pdf")

plt.figure()
im = plt.imshow(np.real(Delta_SE_pseudo_L12), cmap = "RdBu",interpolation="nearest",origin='lower')
#plt.colorbar(label="$\Delta C_E$",fontsize= 22)
plt.xticks(np.arange(0,101,20),fontsize=22) 
plt.yticks(np.arange(0,101,20),fontsize=22)  
# Add a colorbar and set font sizes
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=20)  # Change font size of colorbar ticks
cbar.set_label("$\Delta S_E$", fontsize=22)  # Change font size of the colorbar label

plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/Delta_SE_pseudo_L12.pdf")

plt.figure()
im = plt.imshow(abs(np.real(Delta_SE_pseudo_L12)), cmap = "magma",interpolation="nearest",origin='lower')
#plt.colorbar(label="$\Delta C_E$",fontsize= 22)
plt.xticks(np.arange(0,101,20),fontsize=22) 
plt.yticks(np.arange(0,101,20),fontsize=22)  
# Add a colorbar and set font sizes
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=20)  # Change font size of colorbar ticks
cbar.set_label("$|\Delta S_E|$", fontsize=22)  # Change font size of the colorbar label

plt.show()
plt.savefig("/home/maitri/Documents/krylov_moduler/TFIM_modular/plots/abs_Delta_SE_pseudo_L12.pdf")

'''
