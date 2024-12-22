# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 23:06:47 2024

@author: Utilisateur
"""

import numpy as np
import matplotlib.pyplot as plt 
import random as rnd

def force(r,alpha, beta):
    
    F = (beta/r**8)*((alpha/r**6)-1)*r
    return  F

def LJ(r,epsilon,sigma):
    sigma12 = sigma**12
    sigma6 = sigma**6
    r12= r**12 
    r6 = r**6

    return 4*epsilon*(sigma12/r12 -sigma6/r6)

# epsilon =  34.9
# sigma = 2.78e-10
# alpha = 2*sigma**6
# beta = 24*epsilon*sigma**6

# x = np.linspace(0,1, 10000)
# plt.scatter(x,force(x, beta, alpha))


# epsilon = 34.9*1.3e-23 #Profondeur du puit
epsilon = 100.9*1.3e-23 #Profondeur du puit
sigma= 2.78e-10 #distance sur laquelle le potentiel est effectif
alpha = 2*sigma**6
beta = 24*epsilon*sigma**6


x = np.linspace(5e-11,10e-11, 10000)
pot_LJ = LJ(x, epsilon,sigma)

# print(force_LJ)

plt.plot(x,pot_LJ)


plt.ylim(-10,1e1)
plt.title(f"Potentiel de LJ | $\epsilon$ = {epsilon} ; $\sigma = {sigma}$")
plt.xlabel("Distance r entre deux particules")
plt.ylabel("Potentiel")

plt.show()

#--------------------

plt.plot(x,np.zeros(len(x)), c= "r")

force_LJ = force(x,alpha, beta)
plt.plot(x,force_LJ)

plt.ylim(-1,2)
plt.title(f"Force dérivée de LJ | $\epsilon$ = {epsilon} ; $\sigma = {sigma}$")
plt.xlabel("Distance r entre deux particules")
plt.ylabel("Force")

plt.show()



#g = np.meshgrid(a,b,sparse="True")