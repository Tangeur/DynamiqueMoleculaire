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

df = DataFrames(r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Simulation pertinentes\RSEAUC~1\N_25_T~1.HDF")

df.calcul_temperature(0)
df.calcul_temperature(10)

# df.calcul_energie()

sigma = 2.78e-10
df.calcul_rdf(d_max=3*sigma, nb_points=2000,i_frame=15)

#%%

parameters = {"saving_rate":200,
              "nb_iteration":20_000,
              
              "dimensions":(0.8e-9,0.8e-9),
              "nb_atomes": 207,
              "mass":2e-26, 
              "epsilon": 34.9*1.3e-23,
              "sigma":2.78e-10,
              "dt" : 1e-18,
              "temperature":100,
              "force_Rcut": 2 * 2.78e-10,
              
              "dmin_beetween_atomes":0.6*2.78e-10,
              "grid_type":"Square"}

sigma = 2.78e-10

syst = System(dimensions = (2.1e-9,2.1e-9), nb_atomes = 50, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= sigma, 
              dt = 1e-16, temperature = 100, force_Rcut = 2*sigma)



syst.initialise_system(dmin_beetween_atomes=0.4*2.78e-10,grid_type="Random")
syst.save_initialisation(save_path = path+"/Simulations")

# syst.plot_system()

# syst.add_atome(Atome(np.array([1e-9,1e-9]),syst.defineInitialSpeed(1e-9,1e-9,2e-26,1000,1e-17),2e-26, np.array([0,0])))

syst.plot_system()

syst.iterate(nb_iteration=2000,saving_rate=200,adjust_temperature=False,printing=False)



syst.plot_system()

#---------------------------------------------------------------------

df = DataFrames(syst.path_simulation_file)
df.create_gif(0.0005)
# df.calcul_energie()
df.calcul_temperature(0)
df.calcul_temperature(df.nb_frames-1)

# %%

sigma = 2.78e-10

syst = System(dimensions = (3.4e-9,3.4e-9), nb_atomes = 55, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= sigma, 
              dt = 1e-18, temperature = 200,force_Rcut = 2*sigma)

syst.initialise_system(dmin_beetween_atomes=0.6*2.78e-10,grid_type="Random")
syst.save_initialisation(save_path = path+"/Simulations")

syst.plot_system()

syst.iterate(nb_iteration=600,saving_rate=0,printing=False)
syst.plot_system()

# syst.iterate(nb_iteration=601,saving_rate=200,adjust_temperature=True)

df = DataFrames(syst.path_simulation_file)
df.create_gif(0.0005)
df.calcul_energie()
df.calcul_temperature(df.nb_frames-1)

# df.calcul_rdf(0.2e-7,300)
# df.calcul_rdf(4*2.78e-10,3000)

# %%

#20_000:
# df = DataFrames(r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes\Simulations\N_207_~1.424\N_207_~1.TXT")

#200_000:
df = DataFrames(r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes\Simulations\N_207_~1.063\N_207_~1.TXT")

# df.calcul_rdf(20e-10,300,i_frame=9)
# df.plot_frame(i_frame=9)
# df.create_gif(0.00005)
df.calcul_energie()
df.calcul_temperature(df.nb_frames-1)

#%%


df = DataFrames(r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes\Simulations\N_169_~1.880\N_169_~1.HDF")
df.calcul_rdf(6*sigma, 800,i_frame=0)
# df.plot_frame(0)
