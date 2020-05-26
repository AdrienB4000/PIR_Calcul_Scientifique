from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt

import os
import methodes_resolution as res
import modelisation as mod
import affichage as aff

"""
Todo: bouton méthode

"""

modeles = mod.modeles
noms_modeles = mod.noms_modeles

methodes = res.methodes
noms_methodes = res.noms_methodes

## Paramètres

# Parametres d'affichage
duree = 60 # Duree (jour)
nb_pts_t = 1000
pas_t = duree/nb_pts_t

# Parametres de modelistaion
N_pop = 10000 # Population totale

diffusion = 0 # 0 si non 1 si oui
nb_pts_x = 8


# Choix de u0 pour le cas avec diffusion
u0 = N_pop / (nb_pts_x-1) * np.array([[i,0,(nb_pts_x-1-i),0,0] for i in range(nb_pts_x)])

aff.open(N_pop,duree,pas_t,u0)
