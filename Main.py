# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 21:40:32 2024

@author: Tanguy Ginocchio
"""

import os

#Répertoire dans lequel se trouve le code :
path = r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes"
os.chdir(path)
from Analyse_des_donnees_Dynamique_Moleculaire_2D import DataFrames
from Dynamique_Moléculaire_2D import System


# epsilon =  34.9
# sigma = 2.78e-10
# mass = 2e-26

syst = System(dimensions = (1.8e-9,1.8e-9), nb_atomes = 30, 
              mass = 2e-26, epsilon = 34.9*1.3e-23, sigma= 2.78e-10, 
              dt = 3e-16, temperature = 0)

syst.save_initialisation(save_path = path+"/Simulations")
syst.initialise_system(dmin_beetween_atomes=0.7*2.78e-10,grid="Square")
syst.warning()

syst.plot_system()
syst.calcul_rdf(15e-10,30)

# syst.add_atome(Atome(np.array([20,30]),defineInitialSpeed(20,30,mass,T=100,dt =3e-2),mass, np.array([0,0]),tracking=False,warning=True))

syst.iterate(nb_iteration=2000,saving_rate=200,printing=False)

syst.plot_system()

#---------------------------------------------------------------------

df = DataFrames(syst.path_simulation_file)
df.create_gif()
df.calcul_energie()
df.calcul_rdf(15e-10,30)

##----------------

# df = DataFrames(r"C:\Users\Utilisateur\Desktop\Stockage\Etudes\Université Grenoble Alpes\M2\Etats quantiques de la matière (numérique)\Programmes\Simulations\N_30_T~1.941\N_30_T~1.TXT")
# df.calcul_rdf(15e-10,30)

