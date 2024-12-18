# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 21:40:32 2024

@author: Tanguy Ginocchio
"""

import os
import numpy as np

#Répertoire dans lequel se trouve le code :
path = r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes"
os.chdir(path)
from Analyse_des_donnees_Dynamique_Moleculaire_2D import DataFrames
from Dynamique_Moléculaire_2D import System, Atome



# epsilon =  34.9
# sigma = 2.78e-10
# mass = 2e-26


 #%%
syst = System(dimensions = (1.4e-9,1.4e-9), nb_atomes = 200, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= 2.78e-10, 
              dt = 1e-19, temperature = 100)
syst.save_initialisation(save_path = path+"/Simulations")
# syst.initialise_system(dmin_beetween_atomes=0.15*2.78e-10,grid="Square")
syst.warning()

# syst.add_atome(Atome(np.array([1e-9,1e-9]),syst.defineInitialSpeed(1e-9,1e-9,m =2e-26,T=100,dt =1e-19),2e-26, np.array([0,0]),tracking=False,warning=True))
# syst.add_atome(Atome(np.array([0.1e-9,0.1e-9]),syst.defineInitialSpeed(0.1e-9,0.1e-9,m =2e-26,T=100,dt =1e-19),2e-26, np.array([0,0]),tracking=False,warning=True))
syst.initialise_system(dmin_beetween_atomes=0.15*2.78e-10,grid="Random")

syst.plot_system()
syst.save_frame()

syst.iterate(10)
syst.save_frame()


df = DataFrames(syst.path_simulation_file)
df.calcul_temperature(0)
df.calcul_temperature(1)

#%%
syst = System(dimensions = (0.8e-9,0.8e-9), nb_atomes = 207, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= 2.78e-10, 
              dt = 1e-19, temperature = 0)

syst.save_initialisation(save_path = path+"/Simulations")
syst.initialise_system(dmin_beetween_atomes=0.2*2.78e-10,grid="Square")
syst.warning()

syst.plot_system()

# syst.add_atome(Atome(np.array([20,30]),defineInitialSpeed(20,30,mass,T=100,dt =3e-2),mass, np.array([0,0]),tracking=False,warning=True))

syst.iterate(nb_iteration=8_000,saving_rate=200,printing=False)

syst.plot_system()

#---------------------------------------------------------------------

df = DataFrames(syst.path_simulation_file)
df.create_gif(0.0005)
df.calcul_energie()

# %%

syst = System(dimensions = (1.4e-9,1.4e-9), nb_atomes = 207, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= 2.78e-10, 
              dt = 3e-16, temperature = 200)

syst.save_initialisation(save_path = path+"/Simulations")
syst.initialise_system(dmin_beetween_atomes=0.25*2.78e-10,grid="Random")

syst.plot_system()

syst.iterate(nb_iteration=2_000,saving_rate=200,printing=False)
syst.plot_system()

df = DataFrames(syst.path_simulation_file)
df.create_gif(0.0005)
df.calcul_energie()
df.calcul_temperature(df.nb_frames-1)

# df.calcul_rdf(0.2e-7,300)
# df.calcul_rdf(4*2.78e-10,3000)

# %%

df = DataFrames(r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes\Simulations\N_207_~1.938\N_207_~1.TXT")
# df.calcul_rdf(20e-10,300,i_frame=9)
# df.plot_frame(i_frame=9)
# df.create_gif(0.00005)
df.calcul_energie()
df.calcul_temperature(0)

