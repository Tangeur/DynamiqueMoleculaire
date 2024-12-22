# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 00:34:23 2024

@author: Tanguy GINOCCHIO
"""

import re
import numpy as np
import matplotlib.pyplot as plt 
from decimal import Decimal
import imageio.v2 as imageio
import os
from datetime import datetime
import h5py

kb = 1.38e-23 #J/K

def format_e(n):
    """Format n en écriture scientifique"""
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

class DataFrames :
    """Objet qui charge les données relatives à la simulation dans : 
        self.data_frame[numéro de la simulation][numéro de l'atome][0:positions|1:old_positions][0:x|1:y]"""
    
    def __init__ (self, path_simulation):
        self.path_simulation = path_simulation
        self.path_folder = os.path.dirname(self.path_simulation)
        
        print("Importation des données...")
        
        file = h5py.File(self.path_simulation, "r")
        
        self.dimensions = file.attrs["dimensions"]
        self.mass = file.attrs["mass"]
        self.dt = file.attrs["dt"]
        self.sigma = file.attrs["sigma"]
        self.epsilon = file.attrs["epsilon"]
        self.temperature = file.attrs["temperature"]
        self.force_Rcut = file.attrs["force_Rcut"]
        
        self.simulation_name = file.attrs["simulation_name"]
        
        self.data_frame_positions = []
        self.data_frame_old_positions = []
        
        
        self.nb_frames = len(file["list_frames"])
        self.nb_atomes = file.attrs["nb_atomes"]
        print(f"nb_frames : {self.nb_frames}")
        
        for i in range(self.nb_frames):
            frame_name = file["list_frames"][i].decode("utf-8")
            self.data_frame_positions.append(file[frame_name]["position"][:].tolist())
            self.data_frame_old_positions.append(file[frame_name]["old_position"][:].tolist())
        file.close()
    
        print("Fichier chargé.")

    def create_gif(self,duration = 0.01,title=""):
        """Créé un gif à partir de toutes les frames générée par la simulation"""
        
        if not os.path.exists(self.path_folder+"/all_images"):
            os.mkdir(self.path_folder+"/all_images")

        print("Création des iamges...")

        files_names = []
        images = []
        for i_frame in range(self.nb_frames):
            filename = f"DynamiqueMoléculaire_temp_{i_frame}.png"
            X=np.array([])
            Y=np.array([])
            
            for i_atome in range(self.nb_atomes):
                X = np.append(X,self.data_frame_positions[i_frame][i_atome][0])
                Y = np.append(Y,self.data_frame_positions[i_frame][i_atome][1])
                
            cmap = plt.colormaps['viridis'].resampled(len(X))
            plt.scatter(X,Y,c=range(len(X)),cmap=cmap,s=8)
            
            if title == "":
                plt.title(self.simulation_name)
            else : 
                plt.title(title)
            plt.xlabel("x")
            plt.ylabel("y")
            plt.xlim(0,self.dimensions[0])
            plt.ylim(0,self.dimensions[1])
            
            plt.savefig(self.path_folder + "/all_images/"+filename, dpi = 240)
            
            plt.cla()
            files_names.append(self.path_folder + "/all_images/"+filename)
            images.append(imageio.imread(self.path_folder + "/all_images/"+filename))
            print(f"Gif, frame : {i_frame+1}/{self.nb_frames}")
            
        try :
            print("Compilation du gif...")
            imageio.mimsave(self.path_folder + "/" + self.simulation_name +".gif", images, duration=duration)
            print("Gif exporté sous : "+ self.path_folder +"/" + self.simulation_name +".gif")
        except OSError as error:
            print("ERROR :" + error)
        

    def calcul_energie(self):
        """Créé un graphique représentant l'énergie cinétique dans le système en fonctoin de la frame"""
        
        print("Calcul de l'énegie cinétique au cours de la simulation")
        Ec_list = []
        for i_frame in range(self.nb_frames):
            Ec = 0
            for i_atome in range(self.nb_atomes):
                Ec += 0.5*self.mass*(np.sqrt((self.data_frame_positions[i_frame][i_atome][0]-self.data_frame_old_positions[i_frame][i_atome][0])**2+(self.data_frame_positions[i_frame][i_atome][1]-self.data_frame_old_positions[i_frame][i_atome][1])**2)/self.dt)**2
            Ec_list.append(Ec/self.nb_atomes)
            
        plt.title("Ec au cours des itération")
        plt.ylabel("$E_c$")
        plt.xlabel("Frame")
        plt.plot(range(len(Ec_list)),Ec_list)
        
        date = datetime.now().isoformat().replace(":","-")
        date = date.replace("T", " Time_")
        
        plt.savefig(self.path_folder + f"/Ec au cours des intération Date_{date}.png", dpi = 240)
        plt.show()
        
    def calcul_temperature(self, i_frame):
        
        v2 = 0
        for i_atome in range(self.nb_atomes):
            v2 += (np.sqrt((self.data_frame_positions[i_frame][i_atome][0]-self.data_frame_old_positions[i_frame][i_atome][0])**2+(self.data_frame_positions[i_frame][i_atome][1]-self.data_frame_old_positions[i_frame][i_atome][1])**2)/self.dt)**2
        v2 = v2/self.nb_atomes
        T = v2*self.mass/(2*kb) #On a 1/2 parce qu'on a un système 2D
        print(f"Temperature dans le système : {T}K")
    
    def calcul_rdf(self, d_max, nb_points, i_frame = 0):
        """Calcul la RDF dans le système."""
        
        Lx,Ly = self.dimensions[0], self.dimensions[1]
        dr = d_max/nb_points
        g_r = np.zeros(nb_points)
        d_max_2 = d_max**2
        
        list_atomes = [self.data_frame_positions[i_frame][i] for i in range(self.nb_atomes)]
        
        for atome1 in list_atomes:
            for atome0 in list_atomes:
                dx = atome1[0]-atome0[0]
                dy = atome1[1]-atome0[1]
                
                #print(f"(dx,dy) : ({dx},{dy})")
                
                dx = dx - int(dx*2/Lx)*Lx
                dy = dy - int(dy*2/Ly)*Ly
                r2 = (dx)**2+(dy)**2 #r^2
                
                
                if r2 < d_max_2:
                    r = np.sqrt(r2)
                    g_r[int(r/dr)] += 1
        
        x = np.linspace(0,d_max,nb_points)
        plt.plot(x,g_r)
        plt.show()
        pass
    
    def plot_frame(self, i_frame=0):
        X=np.array([])
        Y=np.array([])
        
        for i_atome in range(self.nb_atomes):
            X = np.append(X,self.data_frame_positions[i_frame][i_atome][0])
            Y = np.append(Y,self.data_frame_positions[i_frame][i_atome][1])
            
        cmap = plt.colormaps['viridis'].resampled(len(X))
        plt.scatter(X,Y,c=range(len(X)),cmap=cmap,s=8)
        
        plt.title(self.simulation_name)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim(0,self.dimensions[0])
        plt.ylim(0,self.dimensions[1])
        
        plt.show()
        # plt.cla()




