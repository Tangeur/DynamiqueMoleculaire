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

#%% Lancement d'une simulation quelconque selon une distribution random

sigma = 2.78e-10

syst = System(dimensions = (2.1e-9,2.1e-9), nb_atomes = 60, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= sigma, 
              dt = 1e-16, temperature = 100, force_Rcut = 2*sigma)



syst.initialise_system(dmin_beetween_atomes=0.4*2.78e-10,grid_type="Random")
syst.save_initialisation(save_path = path+"/Simulations")

syst.plot_system()

syst.iterate(nb_iteration=2000,saving_rate=200,adjust_temperature=True,printing=False)

syst.plot_system()

#-----------------------[Analyse des données]----------------------

df = DataFrames(syst.path_simulation_file)
df.create_gif(0.005)
df.calcul_energie()

df.calcul_temperature(0)
df.calcul_rdf(5*sigma,400,i_frame=0)

df.calcul_temperature(df.nb_frames-1)
df.calcul_rdf(5*sigma,400,i_frame=df.nb_frames-1)

# %% Lancement d'une simulation quelconque selon une distribution carrée

sigma = 2.78e-10

syst = System(dimensions = (1.4e-9,1.4e-9), nb_atomes = 55, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= sigma, 
              dt = 1e-16, temperature = 200,force_Rcut = 2*sigma)

syst.initialise_system(dmin_beetween_atomes=0.6*2.78e-10,grid_type="Square")
syst.save_initialisation(save_path = path+"/Simulations")

syst.plot_system()

syst.iterate(nb_iteration=2_000,saving_rate=200,adjust_temperature=True,printing=False)

syst.plot_system()

#-----------------------[Analyse des données]----------------------

df = DataFrames(syst.path_simulation_file)
df.create_gif(0.005)
df.calcul_energie()

df.calcul_temperature(0)
df.calcul_rdf(5*sigma,400,i_frame=0)

df.calcul_temperature(df.nb_frames-1)
df.calcul_rdf(5*sigma,400,i_frame=df.nb_frames-1)

#%% Reprise d'un fichier pour analyser les données

#On met le chemin du fichier .hdf5
df = DataFrames(r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes\Simulations\N_169_~1.880\N_169_~1.HDF")
df.calcul_rdf(6*sigma, 800,i_frame=0)
df.plot_frame(0)
