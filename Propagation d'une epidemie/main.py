from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt

import os
os.chdir("C:/Users/Anthony/Desktop/COV_A5/")
import methodes_resolution as res
import modelisation as mod
import affichage as aff

modeles = mod.modeles
noms_modeles = mod.noms_modeles

methodes = res.methodes
noms_methodes = res.noms_methodes

## Paramètres

# Parametres d'affichage
duree = 60 # Duree (jour)
nb_pts_t = 100
pas_t = duree/nb_pts_t
temps = np.linspace(0,duree,nb_pts_t)

# Parametres de modelistaion
N_pop = 10000 # Population totale

diffusion = 0 # 0 si non 1 si oui
nb_pts_x = 8

# Methode de resolution
methode = res.Euler_implicite

# Choix de u0 pour le cas avec diffusion
u0 = N_pop / (nb_pts_x-1) * np.array([[i,0,(nb_pts_x-1-i),0] for i in range(nb_pts_x)])

aff.open(methode,N_pop,duree,pas_t,temps,diffusion,u0)








"""
parametres=[N_pop,1.4,0.1,0.2,0.01]
U=methode(u0,mod.f_SIR_D,parametres,temps,pas_t)
dist=np.array(range(nb_pts_x))/(nb_pts_x-1)
for i in range(8):
    plt.plot(dist,U[i][:,0],label=str(i))
plt.legend()
plt.show()#"""

