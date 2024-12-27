# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:09:39 2024

@author: Tanguy GINOCCHIO
"""

import numpy as np
import matplotlib.pyplot as plt 
import random as rnd
from datetime import datetime
import os
import h5py
import sys

rnd.seed(20)
kb = 1.38e-23 #J/K


def format_e(n):
    """Format n en écriture scientifique"""
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

class Atome : 
    
    def __init__(self, position, old_position, mass, force, type = None):
        self.position = position
        self.old_position = old_position
        self.mass = mass
        self.force = force
        self.type = type
    
    def move(self, system, dt) :
        """Fait bouger l'atome en fonction de sa position, son ancienne position et de la résultante des forces qu'il subit des autres atomes."""
        
        Lx, Ly = system.dimensions[0],system.dimensions[1]
        Lx_2, Ly_2 = Lx/2, Ly/2
        abs_Lx_2, abs_Ly_2 = abs(Lx_2), abs(Ly_2)
        
        new_position = 2*self.position - self.old_position + (self.force/self.mass)*dt**2
     
        x_new_relative_position, y_new_relative_position = new_position[0]-Lx_2, new_position[1]-Ly_2
     
        x_n = np.sign(x_new_relative_position)*(abs(x_new_relative_position)//abs_Lx_2)
        y_n = np.sign(y_new_relative_position)*(abs(y_new_relative_position)//abs_Ly_2)
        
        #Tracking :
        #print(f"(x_n,y_n) = ({x_n},{y_n})")
        
        #Warning : 
        if (x_n > 1 or y_n > 1):    
            print(f"[WARNING] : (x_n,y_n) = ({x_n},{y_n}) Les atomes traversent plus d'une boite par itération !")
        
        if abs(x_n) > 0:
            x_new_position_update = (new_position[0]-(np.sign(x_n)*(abs(x_n)+1)*abs_Lx_2)) #On replace le points dans la boite de l'autre côté
        else : 
            x_new_position_update = new_position[0]
            
        if abs(y_n) > 0:    
            y_new_position_update = (new_position[1]-(np.sign(y_n)*(abs(y_n)+1)*abs_Ly_2))
        else :
            y_new_position_update = new_position[1]
        
        x_old_position_update = x_new_position_update + - (new_position[0]-self.position[0])
        y_old_position_update = y_new_position_update + - (new_position[1]-self.position[1])
            
        #Tracking :
        #print(f"position : ({self.position[0]},{self.position[1]})\n"
        #      f"old_position : ({self.old_position[0]},{self.old_position[1]})\n"
        #      f"new_position : ({new_position[0]},{new_position[1]})\n"
        #      f"new_position_update : ({x_new_position_update},{y_new_position_update})\n"
        #      f"old_position_update : ({x_old_position_update},{y_old_position_update})")
        
        self.position,self.old_position = np.array([x_new_position_update,y_new_position_update]), np.array([x_old_position_update,y_old_position_update])

class System:
    
    def __init__(self, dimensions, nb_atomes, mass, epsilon, sigma, dt, temperature, force_Rcut):
        self.dimensions = dimensions
        self.nb_atomes = nb_atomes
        self.mass = mass
        self.epsilon = epsilon
        self.sigma = sigma
        self.dt = dt 
        self.temperature = temperature
        self.force_Rcut2 = force_Rcut*force_Rcut
        
        self.alpha = 2*sigma**6
        self.beta = 24*epsilon*sigma**6
        self.iteration = 0
        self.list_atomes = []
        self.list_frames = []
        self.grid_type = None
        
        self.save_path = ""     #Chemin du dossier qui contient les dossiers de toutes les simulations
        self.path_simulation_folder = ""  #Chemin du dossier qui contient les fichiers de la simulation qu'on crée
        self.simulatoin_name = ""   #Nom de la simulation (relatif à ses paramètres)
        self.path_simulation_file = ""  #Chemin du fichier texte qui contient les données de la simulation (position des atomes, etc...)
            
    
    def defineInitialSpeed(self, x_pos, y_pos, m, T, dt):
        """Calcul la old_position que doit avoir un atome, pour avoir une vitesse thermique de température T."""
        
        v_norm = np.sqrt(2*kb*T/m)
        angle = rnd.uniform(0, 2*np.pi)
        
        v = np.array([v_norm*np.cos(angle),v_norm*np.sin(angle)])
        old_position_calculated = np.array([x_pos-v[0]*dt,y_pos-v[1]*dt])
        
        return old_position_calculated
    
    def initialise_system(self,  dmin_beetween_atomes, grid_type = "Random"):
        """Crée le system selon le réseau souhaité selon : Random/Square avec un paramètre de maille : dmin_beetween_atomes."""
        
        print(f"Création d'un système de dimensions {self.dimensions} consitué de {self.nb_atomes} atomes")
        
        if grid_type == "Square" :
            
            print(f"[Info] Initialisation d'un système selon un réseau carré de paramètre de maille : {dmin_beetween_atomes} \nLe nombre d'atomes va être adapté au réseau.")
            
            self.grid_type = "Square"
            
            Lx,Ly = self.dimensions[0],self.dimensions[1]
            
            if Lx%dmin_beetween_atomes != 0:
                Lx = (Lx//dmin_beetween_atomes+1)*dmin_beetween_atomes
            if Ly%dmin_beetween_atomes != 0:
                Ly = (Ly//dmin_beetween_atomes+1)*dmin_beetween_atomes
            # print(f"(Lx,Ly)/dmin_beetween_atomes : ({Lx/dmin_beetween_atomes},{Ly/dmin_beetween_atomes})")
            
            self.dimensions = (Lx,Ly)
            print(f"[Info] Ajustement de la taille du système : (Lx,Ly) = ({Lx},{Ly})")
            
            N_xgrid = int(Lx / dmin_beetween_atomes)
            N_ygrid = int(Ly / dmin_beetween_atomes)
            
            boxes = []
            
            for i in range(N_xgrid):
                for j in range(N_ygrid):
                    boxes.append((i,j))
            
            self.nb_atomes = len(boxes)
            print(f"[Info] Ajustement du nombre d'atomes : {self.nb_atomes}")
            
            for boxe in boxes:
                x_boxe, y_boxe = boxe[0],boxe[1]
                x_pos = x_boxe*dmin_beetween_atomes + dmin_beetween_atomes/2
                y_pos = y_boxe*dmin_beetween_atomes + dmin_beetween_atomes/2
                position = np.array([x_pos,y_pos])
                old_position = self.defineInitialSpeed(x_pos,y_pos,self.mass,self.temperature, self.dt)
                force = np.array((0,0))
                self.list_atomes.append(Atome(position,old_position,self.mass,force))

        
        elif grid_type == "Hexagon":
            print("[ERREUR] Le système ne peut pas être initialisé comme cela.")
            pass
        
        elif grid_type == "Random":
            
            print(f"[Info] Initialisation d'un système dans une distribution aléatoire, avec une distance minimum entre chaque atome de : {dmin_beetween_atomes}")
            
            self.grid_type = "Random"
            
            Lx,Ly = self.dimensions[0],self.dimensions[1]
            
            #On adapte la taille du système pour que les côtés soient un multiple de 'dmin_beetween_atomes'
            if Lx%dmin_beetween_atomes != 0:
                Lx = (Lx//dmin_beetween_atomes+1)*dmin_beetween_atomes
            if Ly%dmin_beetween_atomes != 0:
                Ly = (Ly//dmin_beetween_atomes+1)*dmin_beetween_atomes
            # print(f"(Lx,Ly)/dmin_beetween_atomes : ({Lx/dmin_beetween_atomes},{Ly/dmin_beetween_atomes})")
            
            self.dimensions = (Lx,Ly)
            print(f"Ajustement de la taille du système : (Lx,Ly) = ({Lx},{Ly})")
            
            N_xgrid = int(Lx / dmin_beetween_atomes)
            N_ygrid = int(Ly / dmin_beetween_atomes)
            
            boxes = []
            
            for i in range(N_xgrid):
                for j in range(N_ygrid):
                    boxes.append((i,j))
                    
            #print(f"Vérif : (N_xgrid,N_ygrid) = ({N_xgrid},{N_ygrid}) | N_xgrid*N_ygrid = {N_xgrid*N_ygrid} | len(boxes) = {len(boxes)}" )
            
            for i in range(N_xgrid):
                boxes.remove((i,N_ygrid-1))
                
            for j in range(N_ygrid-1):
                boxes.remove((N_xgrid-1,j))
            
            for i in range(self.nb_atomes):
                
                if len(boxes) == 0:
                    print(f"On a pu ajouter {len(self.list_atomes)} atomes.")
                    print("[ERROR] : On ne peut pas rajouter d'atomes sur la grille... Veuillez adapter le nombre d'atomes dans le système.")
                    raise ValueError(1) #On arrête le programme car on n'a pas pu générer toute la grille...
                    
                
                boxe = boxes[rnd.randint(0,len(boxes)-1)]
                x_boxe, y_boxe = boxe[0], boxe[1]
                
                x_pos = x_boxe*dmin_beetween_atomes + rnd.uniform(0, dmin_beetween_atomes)
                y_pos = y_boxe*dmin_beetween_atomes + rnd.uniform(0, dmin_beetween_atomes)
                position = np.array([x_pos,y_pos])
                old_position = self.defineInitialSpeed(x_pos,y_pos,self.mass,self.temperature, self.dt)
                force = np.array((0,0))
                self.list_atomes.append(Atome(position,old_position,self.mass,force))
                
                for i in range(-1,2,1):
                    for j in range(-1,2,1):
                        if (i+x_boxe,j+y_boxe) in boxes:
                            boxes.remove((i+x_boxe,j+y_boxe))   
                
        else : 
            print("[ERREUR] Impossible de créer le réseau")
              
    def force(self, atome0,atome1,alpha,beta,Lx,Ly,force_Rcut2):
        """Calcul la force exercé par l'atome0 sur l'atome 1."""
        
        dx = atome1.position[0]-atome0.position[0]
        dy = atome1.position[1]-atome0.position[1]
        
        #print(f"(dx,dy) : ({dx},{dy})")
        
        dx = dx - int(dx*2/Lx)*Lx   #Application des conditions aux bords périodiques
        dy = dy - int(dy*2/Ly)*Ly
        r2 = (dx)**2+(dy)**2 #r^2
        
        if r2 < self.force_Rcut2 : #Si l'atome est au délà du force_Rcut on ne le prend pas en compte 
            if (r2==0):
                return np.array([.0,.0])
            
            r6 = r2**3 #r^6
            r8 = r2**4 #r^8
            # r = np.sqrt((dx)**2+(dy)**2)
            beta_r8 =beta/r8
            alpha_r6 = alpha/r6 - 1
            
            F_x = beta_r8*alpha_r6*dx
            F_y = beta_r8*alpha_r6*dy
            # F = (beta/r**8)*((alpha/r**6)-1)*(atome1.position-atome0.position)
            return  np.array([F_x,F_y])
        else :
            return np.array([.0,.0])
    
    def calculate_forces(self):
        """Calcul la force que subit chaque atome du à la présence de tous les autres atomes."""
        
        for atome1 in self.list_atomes:
            force_resultante = np.array([.0,.0])
            for atome0 in self.list_atomes:
                force_resultante += self.force(atome0,atome1,self.alpha,self.beta,self.dimensions[0],self.dimensions[1],self.force_Rcut2)
            atome1.force = force_resultante
    
    def plot_system(self):
        """Crée et affiche un graph des particules dans le système."""
        
        X=np.array([])
        Y=np.array([])
        
        for atome in self.list_atomes:
            X = np.append(X,atome.position[0])
            Y = np.append(Y,atome.position[1])
        plt.xlim(0,self.dimensions[0])
        plt.ylim(0,self.dimensions[1])
        
        cmap = plt.colormaps['viridis'].resampled(len(X))
        plt.scatter(X,Y,c=range(len(X)),cmap=cmap,s=8)
        plt.show()
    
    def iterate(self, nb_iteration, saving_rate = 0, adjust_temperature = False, printing=False):
        """Fait passer le système au pas de temps suivant.
        
        -saving_rate = 10 :  Sauvegarde les données dans le fichier toutes les 10 frames
        -adjust_temperature = True : Corrige la vitesse des particules pour qu'elle reste à la vitesse thermique lié à la température du système
        -printing = True :  Affiche les plot au fur et à mesure. """

        print("Itération du système...")
        for i in range(nb_iteration):
            
            self.calculate_forces() #Calcul des forces qui s'exercent sur chaque atome 
            for atome in self.list_atomes:
                atome.move(self, self.dt) #On fait bouger les atomes d'un pas de temps
            
            if adjust_temperature == True: #On ajuste la vitesse des atomes pour qu'ils gardent la même vitesse thermique 
                thermal_speed_norm = np.sqrt(2*kb*self.temperature/self.mass) #norme de la vitesse thermique
                for atome in self.list_atomes:
                    vect =  atome.old_position - atome.position
                    vect_norm = np.sqrt(vect[0]**2+vect[1]**2)  #Norme du vecteur allant de la position vers la old_position                  
                    a = thermal_speed_norm*self.dt/vect_norm    #Facteur de proportionalité 
                    new_old_position = atome.position + vect * a 
                    atome.old_position = new_old_position

            self.iteration += 1
            
            try : 
                if i % saving_rate == 0: 
                    if printing:
                        self.plot_system() 
                    
                    self.save_frame() #On sauvegarde la simulation 
                    print(f"Implementation {i}/{nb_iteration}")
                    
            except ZeroDivisionError :
                #On ne fait rien, on ne sauvegarde pas de trace de l'évolution du système
                pass
        
        file = h5py.File(self.path_simulation_folder+"/"+self.simulation_name+".hdf5",'a')        
        file.create_dataset(name="list_frames", data=self.list_frames) #On enregistre la liste des frames dans le fichier pour la lecture des données 
        file.close()
        print("Séquence d'itération terminée")
    
    def add_atome(self, Atome):
        """Ajoute un atome au système."""
        
        self.list_atomes.append(Atome)
    
    def save_frame(self): 
        """Sauvegarde les données de la frame dans le fichier de simulation"""
        
        positions = []
        old_positions = []
        
        self.list_frames.append(f"frame_{self.iteration}")
        
        for atome in self.list_atomes:
            positions.append(atome.position.tolist())
            old_positions.append(atome.old_position.tolist())     
        
        file = h5py.File(self.path_simulation_folder+"/"+self.simulation_name+".hdf5",'a')        
        
        group = file.create_group(f"frame_{self.iteration}")
        file.flush()
        group.create_dataset("position",data=positions)
        group.create_dataset("old_position",data=old_positions)
        group.attrs.create("iteration", self.iteration)
        file.close()
    
    def save_initialisation(self, save_path, simulation_name = ""):
        
        """Cette fonction  : 
            ->Vérifie que le système a bien été initialisé 
            ->Crée le dossier dans lequel on enregistre
            ->Crée le fichier qui contient les données que la simulation va enregistrer
            
            Le nom de la fonction se génère tout seul en fonction des paramètres du système.
            """
        
        if self.grid_type == None:
            print("[ERREUR] Vous n'avez pas initialisé le système avec : 'initialise_system()'")
            sys.exit()
        
        self.list_frames = []
        
        #Assure du dossier dans lequel sont enregistré tous les dossier de simulation 
        self.save_path = save_path
        if not os.path.isdir(self.save_path):
            os.mkdir(self.save_path)
        
        date = datetime.now().isoformat().replace(":","-")
        date = date.replace("T", " Time_")
        date = date.split(".")[0]
        
        #Définit le nom de la simulation
        if simulation_name == "":
            self.simulation_name = f"N={self.nb_atomes} T={self.temperature}K epsilon={format(self.epsilon,'.4g')}kBJ sigma={format(self.sigma,'.4g')}pm dimensions=({format(self.dimensions[0],'.4g')},{format(self.dimensions[1],'.4g')}) dt={format(self.dt,'.4g')}" +f" Date_{date}"
        else :
            self.simulation_name = simulation_name
        
        #Création du dossier de la simulation
        self.path_simulation_folder = self.save_path+"/"+self.simulation_name
        os.mkdir(self.path_simulation_folder)
        
        self.path_simulation_file = self.path_simulation_folder+"/"+self.simulation_name+".hdf5"
        
        #Création et ouverture du fichier de simulation
        file = h5py.File(self.path_simulation_file, 'a')

        #Enregistre les attributs de la simulation dans le fichier
        file.attrs.create("dimensions", self.dimensions)
        file.attrs.create("nb_atomes",self.nb_atomes)
        file.attrs.create("mass",self.mass)
        file.attrs.create("epsilon", self.epsilon)
        file.attrs.create("sigma",self.sigma)
        file.attrs.create("dt", self.dt)
        file.attrs.create("temperature", self.temperature)
        file.attrs.create("force_Rcut",np.sqrt(self.force_Rcut2))
        file.attrs.create("simulation_name",self.simulation_name)

        file.flush()
        file.close()





