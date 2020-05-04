from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

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
duree = 80 # Duree (jour)
nb_pts = 10000
pas = duree/nb_pts
temps = np.linspace(0,duree,nb_pts)

# Parametres de modelistaion
N = 100000 # Population totale
beta_0 = 1.4 # Contacts par pers. infectee par jour
gamma_0 = 1/12 # Inverse de la duree de guerison
Lambda = 0.012/365*N # Natalite (nb/jour)
mu = 0.009/365 # Taux de mortalite (nb/pers/jour)
alpha_0 = 1/5 # Inverse de la duree d'incubation
# R_0 represente le nombre de contact par pers. infectee, varie selon le choix de f.
D_0 = 1

parametres = [N,beta_0,gamma_0,Lambda,mu,alpha_0,D_0]

# Parametres de resolution
E_0 = 0
I_0 = 10
S_0 = N-I_0
R_0 = 0
u_0 = np.array([S_0,E_0,I_0,R_0]) # [S,E,I,R]

methode = res.Runge_Kutta_4 # Methode de resolution
modele = mod.f_SIR # Modele d'epidemie
cree = True # Affiche le resultat sur une nouvelle fenetre

aff.open(methode,N,duree,pas,temps,Lambda,mu,D_0)