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
from decimal import Decimal

rnd.seed(20)
kb = 1.38e-23 #J/K


def format_e(n):
    """Format n en écriture scientifique"""
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

class Atome : 
    
    def __init__(self, position, old_position, mass, force, type = None, tracking = False, warning = False):
        self.position = position
        self.old_position = old_position
        self.mass = mass
        self.force = force
        self.type = type
        self.warning = warning
    
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
        
        if self.warning:
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
    
    def __init__(self, dimensions, nb_atomes, mass, epsilon, sigma, dt, temperature):
        self.dimensions = dimensions
        self.nb_atomes = nb_atomes
        self.mass = mass
        self.epsilon = epsilon
        self.sigma = sigma
        self.dt = dt 
        self.temperature = temperature
        
        self.alpha = 2*sigma**6
        self.beta = 24*epsilon*sigma**6
        self.iteration = 0
        self.list_atomes = []
        
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
    
    def initialise_system(self,  dmin_beetween_atomes, grid = "Random"):
        """Créé le system selon le réseau souhaité selon : Random/Square/Hexagon avec un paramètre de maille : dmin_beetween_atomes."""
        
        print(f"Création d'un système de dimensions {self.dimensions} consitué de {self.nb_atomes} atomes")
        
        if grid == "Square" :
            
            Lx,Ly = self.dimensions[0],self.dimensions[1]
            
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
            
            self.nb_atomes = len(boxes)
            print(f"Ajustement du nombre d'atomes : {self.nb_atomes}")
            
            for boxe in boxes:
                x_boxe, y_boxe = boxe[0],boxe[1]
                x_pos = x_boxe*dmin_beetween_atomes + dmin_beetween_atomes/2
                y_pos = y_boxe*dmin_beetween_atomes + dmin_beetween_atomes/2
                position = np.array([x_pos,y_pos])
                old_position = self.defineInitialSpeed(x_pos,y_pos,self.mass,self.temperature, self.dt)
                force = np.array((0,0))
                self.list_atomes.append(Atome(position,old_position,self.mass,force))

        
        elif grid == "Hexagon":
            print("[ERREUR] Le système ne peut pas être initialisé comme cela.")
            pass
        
        elif grid == "Random":
            
            Lx,Ly = self.dimensions[0],self.dimensions[1]
            
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
                    print("[ERROR] : On ne peut pas rajouter d'atomes sur la grille...")
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
              
    def force(self, atome0,atome1,alpha,beta,Lx,Ly):
        """Calcul la force exercé par l'atome0 sur l'atome 1."""
        
        dx = atome1.position[0]-atome0.position[0]
        dy = atome1.position[1]-atome0.position[1]
        
        #print(f"(dx,dy) : ({dx},{dy})")
        
        dx = dx - int(dx*2/Lx)*Lx
        dy = dy - int(dy*2/Ly)*Ly
        r2 = (dx)**2+(dy)**2 #r^2
        
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
    
    def calculate_forces(self):
        """Calcul la force que subit chaque atome du à la présence de tous les autres atomes."""
        
        for atome1 in self.list_atomes:
            force_resultante = np.array([.0,.0])
            for atome0 in self.list_atomes:
                force_resultante += self.force(atome0,atome1,self.alpha,self.beta,self.dimensions[0],self.dimensions[1])
            atome1.force = force_resultante
    
    def plot_system(self):
        """Crée un graph des particules du système."""
        
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
      
    def save_frame(self): 
        folder = open(self.path_simulation_folder+"/"+self.simulation_name+".txt", "a")
        folder.write(f"#Itération : {self.iteration}\n")
        for atome in self.list_atomes:
            folder.write(str(atome.position)+","+str(atome.old_position)+"\n")
        folder.close()
    
    def iterate(self, nb_iteration, saving_rate = 0, printing=False):
        """Fait passer le système au pas de temps suivant.
        printing : (True) Affiche les plot au fur et à mesure. """

        for i in range(nb_iteration):
            self.calculate_forces()
            for atome in self.list_atomes:
                atome.move(self, self.dt)
            self.iteration += 1
            
            try : 
                if i % saving_rate == 0:
                    if printing:
                        self.plot_system() 
                    
                    self.save_frame()
                    print(f"Implémentation {i}/{nb_iteration}")
                    
            except ZeroDivisionError :
                #On ne fait rien, on ne sauvegarde pas de trace de l'évolution du système
                pass
    
    
    
    def add_atome(self, Atome):
        """Ajoute un atome au système."""
        
        self.list_atomes.append(Atome)
    
    def warning(self):
        """Active l'affichage des problèmes rencontrés lors de la simulation (Ex: lorsqu'une particule traverse plus d'une boite)"""
        for i in self.list_atomes:
            i.warning = True
    
    def save_initialisation(self, save_path, simulation_name = ""):
        self.save_path = save_path
        if not os.path.isdir(self.save_path):
            os.mkdir(self.save_path)
        
        date = datetime.now().isoformat().replace(":","-")
        date = date.replace("T", " Time_")
        
        if simulation_name == "":
            self.simulation_name = f"N={self.nb_atomes};T={self.temperature}K;epsilon={self.epsilon}kBJ;sigma={self.sigma}pm;dimensions=({self.dimensions[0]},{self.dimensions[1]});dt={self.dt}" +f" Date_{date}"
        else :
            self.simulation_name = simulation_name
        
        self.path_simulation_folder = self.save_path+"/"+self.simulation_name
        os.mkdir(self.path_simulation_folder)
        
        self.path_simulation_file = self.path_simulation_folder+"/"+self.simulation_name+".txt"
        
        folder = open(self.path_simulation_folder+"/"+self.simulation_name+".txt", "a")
        folder.write("#Position : [x,y],Old_position:[old_x,old_y]\n")
        folder.close()




